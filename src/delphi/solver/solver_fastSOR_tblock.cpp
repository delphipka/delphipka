// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: new file, temporal blocking of the SOR half sweeps.

/*
 * solver_fastSOR_tblock.cpp
 *
 * Several red-black half sweeps in one pass over the grid. The plain loop of itit() streams phimap1, phimap2 and the salt map
 * from memory once per half sweep; here a half sweep s runs on plane n as soon as the sweep s-1 has finished plane n+1, so
 * the planes a sweep needs are still in the cache when the next sweep arrives and the arrays cross the memory bus once per
 * block of sweeps instead of once per sweep. Every point is updated by exactly the same operations in the same order as in
 * itrOddPoints()/itrEvenPoints() (stencil, dielectric-boundary term, face reset, charge term), so the result is bit-identical
 * to the plain loop. The convergence test of itit() needs the grid after two particular iterations, so those are left to the
 * plain loop (see tbIterations()). Used for one thread per process only (the MPI ranks), without periodic boundaries, with
 * a fixed relaxation factor and without the nonlinear terms; otherwise itit() runs the plain loop.
 * DELPHI_TB=0 switches it off, DELPHI_TB_DEPTH sets the half sweeps per block (default 8).
 */

#include "solver_fastSOR.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

namespace
{
    // one plane of a half sweep: the three variants of the loops in itrOddPoints()/itrEvenPoints(), same operation order
    // dst = the array being updated, src = the other colour, a = first - stencil offsets of the colour (odd: -1 / even: +1)
    template <int MODE>
    inline void planeSweep(delphi_real* __restrict dst, const delphi_real* __restrict src, const delphi_real* __restrict sf,
                           delphi_integer star, delphi_integer fin, bool odd, delphi_integer lat1, delphi_integer lat2,
                           delphi_integer long1, delphi_integer long2, delphi_real om1, delphi_real sixth)
    {
        if (odd)
        {
            for (delphi_integer ix = star; ix <= fin; ix++)
            {
                const delphi_real temp1 = src[ix-1]         + src[(ix-1)-1];
                const delphi_real temp2 = src[(ix-1)+lat1]  + src[(ix-1)-lat2];
                const delphi_real temp3 = src[(ix-1)+long1] + src[(ix-1)-long2];
                if (MODE == 0) dst[ix-1] = dst[ix-1]*om1 + (temp1+temp2+temp3)*sf[ix-1];
                else           dst[ix-1] = dst[ix-1]*om1 + (temp1+temp2+temp3)*sixth;
            }
        }
        else
        {
            for (delphi_integer ix = star; ix <= fin; ix++)
            {
                const delphi_real temp1 = src[ix-1]         + src[(ix-1)+1];
                const delphi_real temp2 = src[(ix-1)+lat2]  + src[(ix-1)-lat1];
                const delphi_real temp3 = src[(ix-1)+long2] + src[(ix-1)-long1];
                if (MODE == 0) dst[ix-1] = dst[ix-1]*om1 + (temp1+temp2+temp3)*sf[ix-1];
                else           dst[ix-1] = dst[ix-1]*om1 + (temp1+temp2+temp3)*sixth;
            }
        }
    }
}

// the plane (1..iGrid-2) whose updated range [sta,fin] holds the 1-based position ix, 0 if none
static int planeOf(const vector<delphi_integer>& sta, const vector<delphi_integer>& fin, delphi_integer iGrid, delphi_integer ix)
{
    // fin[] grows with the plane number
    int lo = 1, hi = (int)iGrid - 2;
    while (lo < hi)
    {
        const int mid = (lo + hi)/2;
        if (fin[mid] < ix) lo = mid + 1; else hi = mid;
    }
    return (lo >= 1 && lo <= (int)iGrid - 2 && sta[lo] <= ix && ix <= fin[lo]) ? lo : 0;
}

// entries[e] (a 1-based grid position) bucketed by plane, stable; start has iGrid+1 offsets, bucket 0 is "in no plane"
static void bucketByPlane(const vector<delphi_integer>& pos, const vector<delphi_integer>& sta, const vector<delphi_integer>& fin,
                          delphi_integer iGrid, vector<int>& start, vector<int>& order, bool& inPlane)
{
    const int G = (int)iGrid;
    start.assign(G + 2, 0);
    vector<int> pl(pos.size());
    for (size_t e = 0; e < pos.size(); e++)
    {
        pl[e] = planeOf(sta, fin, iGrid, pos[e]);
        if (0 == pl[e]) inPlane = false;
        start[pl[e] + 1]++;
    }
    for (int n = 0; n <= G; n++) start[n+1] += start[n];
    order.assign(pos.size(), 0);
    vector<int> fill(start.begin(), start.end() - 1);
    for (size_t e = 0; e < pos.size(); e++) order[fill[pl[e]]++] = (int)e;
}

bool CDelphiFastSOR::tbPrepare()
{
    tbReady = false;
    const char* env = std::getenv("DELPHI_TB");
    if (env && '0' == env[0]) return false;
    if (!bQmapZero || bFixedRelaxParam || 0.0 < fGridConverge || bLogGraph) return false;
    if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]) return false;
    if (iIterateInterval < 4 || iGrid < 8) return false;
#ifdef PARALLEL_OMP
    if (omp_get_max_threads() != 1) return false;
#endif
    if (phimap1.size() < (size_t)iHalfGridNum || phimap2.size() < (size_t)iHalfGridNum) return false;
    if (bndx1.empty() || bndx3.empty()) return false;

    tbDepth = 8;
    if (std::getenv("DELPHI_TB_DEPTH")) tbDepth = std::max(2, std::atoi(std::getenv("DELPHI_TB_DEPTH")));
    tbDepth &= ~1;     // whole iterations

    const int G = (int)iGrid;
    bool ok = true;

    // dielectric-boundary and charge entries: positions are 1-based, every entry must lie in the updated range of a plane
    vector<delphi_integer> pos;
    for (int c = 0; c < 2; c++)
    {
        const vector<delphi_integer>& sta = c ? sta2 : sta1;
        const vector<delphi_integer>& fin = c ? fi2 : fi1;

        pos.clear();
        for (delphi_integer n = c ? iDielecBndyEven : 0; n < (c ? iDielecBndyOdd : iDielecBndyEven); n++) pos.push_back(prgiBndyDielecIndex[n]);
        bool inPlane = true;
        bucketByPlane(pos, sta, fin, iGrid, tbStart[c][0], tbOrder[c][0], inPlane);
        ok = ok && inPlane && 0 == tbStart[c][0][1];     // bucket 0 must be empty

        pos.clear();
        for (delphi_integer n = c ? iCrgedGridEven : 0; n < (c ? iCrgedGridSum : iCrgedGridEven); n++) pos.push_back(prgiCrgPose[n]);
        inPlane = true;
        bucketByPlane(pos, sta, fin, iGrid, tbStart[c][2], tbOrder[c][2], inPlane);
        ok = ok && inPlane && 0 == tbStart[c][2][1];

        // face reset: entry e = 2*n + w, target position (1-based) iy (w = 0) or iy+(G+1)/2-1 (w = 1) of the serial loops
        pos.clear();
        const delphi_integer star = c ? (G+2)/2 : (G+1)/2, fin2 = c ? (G*(G-1)-1)/2 : (G*(G-1)-2)/2;
        delphi_integer iy = c ? G*(G+2)/2-G+1 : G*(G+1)/2-G+1;
        for (delphi_integer n = 0; n < fin2-star+1; n++)
        {
            iy += G;
            pos.push_back(iy);
            pos.push_back(iy + ((G+1)/2 - 1));
        }
        inPlane = true;
        bucketByPlane(pos, sta, fin, iGrid, tbStart[c][1], tbOrder[c][1], inPlane);   // bucket 0 = outside the updated ranges
        tbResetPos[c] = pos;
    }
    if (!ok) return false;
    tbReady = true;
    return true;
}

// the number of iterations, starting at iteration itr (1-based), that can be run blocked: the ones between the two that the
// convergence test needs (the copy of phimap2 after iteration k*icon1-1 and the comparison after iteration k*icon1)
int CDelphiFastSOR::tbIterations(int itr, int maxItr) const
{
    if (!tbReady) return 0;
    const int pos = itr % iIterateInterval;
    if (pos < 1 || pos > iIterateInterval - 2) return 0;
    int b = iIterateInterval - 1 - pos;
    b = std::min(b, maxItr - itr);        // the plain loop's own limits stay in force near the end
    return b >= 2 ? b : 0;
}

void CDelphiFastSOR::tbHalfSweeps(int count)
{
    const int G = (int)iGrid;
    const bool salt = fZero < std::abs(fIonStrength);
    delphi_real* ph[2]        = { phimap1.data(), phimap2.data() };
    const delphi_real* sf[2]  = { prgfSaltMap1.data(), prgfSaltMap2.data() };
    const vector<delphi_integer>* sta[2] = { &sta1, &sta2 };
    const vector<delphi_integer>* fin[2] = { &fi1, &fi2 };
    const vector<delphi_real>* bnd[2][2]  = { { &bndx1, &bndx2 }, { &bndx3, &bndx4 } };

    for (int base = 0; base < count; base += tbDepth)
    {
        const int t = std::min(tbDepth, count - base);
        for (int p = 1; p < G - 2 + t; p++)
        {
            for (int s = 0; s < t; s++)
            {
                const int n = p - s;
                if (n < 1 || n > G - 2) continue;
                const int c = (base + s) & 1;            // 0: odd points (phimap1), 1: even points (phimap2)
                delphi_real* dst = ph[c];
                const delphi_real* src = ph[1 - c];

                if (salt) planeSweep<0>(dst, src, sf[c], (*sta[c])[n], (*fin[c])[n], 0 == c, lat1, lat2, long1, long2, om1, sixth);
                else      planeSweep<1>(dst, src, sf[c], (*sta[c])[n], (*fin[c])[n], 0 == c, lat1, lat2, long1, long2, om1, sixth);

                // dielectric-boundary terms of this plane
                for (int k = tbStart[c][0][n]; k < tbStart[c][0][n+1]; k++)
                {
                    const delphi_integer e = tbOrder[c][0][k] + (c ? iDielecBndyEven : 0);
                    const delphi_integer ix = prgiBndyDielecIndex[e];
                    const vector<delphi_real>& db = prgfBndyDielec[e];
                    delphi_real temp1, temp2, temp3;
                    if (0 == c)
                    {
                        temp1 = src[(ix-1)-1]*db[0]     + src[ix-1]*db[1];
                        temp2 = src[(ix-1)-lat2]*db[2]  + src[(ix-1)+lat1]*db[3];
                        temp3 = src[(ix-1)-long2]*db[4] + src[(ix-1)+long1]*db[5];
                    }
                    else
                    {
                        temp1 = src[ix-1]*db[0]         + src[ix]*db[1];
                        temp2 = src[(ix-1)-lat1]*db[2]  + src[(ix-1)+lat2]*db[3];
                        temp3 = src[(ix-1)-long1]*db[4] + src[(ix-1)+long2]*db[5];
                    }
                    dst[ix-1] += temp1 + temp2 + temp3;
                }

                // reset the face points the stencil loop overwrote
                for (int k = tbStart[c][1][n]; k < tbStart[c][1][n+1]; k++)
                {
                    const int e = tbOrder[c][1][k];
                    dst[tbResetPos[c][e]-1] = (*bnd[c][e & 1])[e >> 1];
                }

                // charge terms
                for (int k = tbStart[c][2][n]; k < tbStart[c][2][n+1]; k++)
                {
                    const delphi_integer e = tbOrder[c][2][k] + (c ? iCrgedGridEven : 0);
                    dst[prgiCrgPose[e]-1] += prgfCrgValA[e];
                }
            }
        }
    }

    // resets whose target no stencil loop touches (outside every updated range): a no-op, done anyway to match the plain loop
    for (int c = 0; c < 2; c++)
        for (int k = tbStart[c][1][0]; k < tbStart[c][1][1]; k++)
        {
            const int e = tbOrder[c][1][k];
            ph[c][tbResetPos[c][e]-1] = (*bnd[c][e & 1])[e >> 1];
        }
}
