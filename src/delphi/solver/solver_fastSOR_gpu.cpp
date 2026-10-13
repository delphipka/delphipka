// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: new file, optional OpenCL (GPU) version of the iteration loop of itit().

/*
 * solver_fastSOR_gpu.cpp
 *
 * The whole iteration loop of CDelphiFastSOR::itit() on a GPU through OpenCL (gpu_cl.h). The data stays on the device between
 * iterations: one iteration is the odd half sweep and the even half sweep, each made of the kernels that reproduce the loops of
 * itrOddPoints()/itrEvenPoints() in the same order (plane stencil, dielectric-boundary correction, face reset, charge terms), and
 * every iIterateInterval iterations the change of phimap2 is reduced on the device and only two numbers come back to the host.
 * Nothing is written to the host arrays until the loop has finished, so any OpenCL failure just returns false and itit() runs
 * on the CPU as before. Cases the kernels do not cover (periodic boundaries, the nonlinear solver, Chebyshev relaxation,
 * a grid-energy stop, graph logging, a strided convergence test) also return false.
 */

#include "solver_fastSOR.h"
#include "gpu_cl.h"
#include <algorithm>
#include <chrono>
#include <cmath>

bool CDelphiFastSOR::gpuItit(delphi_real maxres, int& itrOut)
{
    using namespace delphicl;

    const char* gpuEnv = std::getenv("DELPHI_GPU");
    if (gpuEnv && '0' == gpuEnv[0]) return false;
    if (!bQmapZero || bFixedRelaxParam || 0.0 < fGridConverge || bLogGraph || 1 != iConvergeFract) return false;
    if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]) return false;

    int minGrid = GPU_MIN_GRID;
    if (std::getenv("DELPHI_GPU_MIN_GRID")) minGrid = std::atoi(std::getenv("DELPHI_GPU_MIN_GRID"));
    else if (gpuEnv && '1' == gpuEnv[0]) minGrid = 0;
    if (iGrid < minGrid) return false;

    Runtime& R = Runtime::get(sizeof(delphi_real) == 8);
    if (!R.ok()) return false;
    const Api& A = R.api;

    const int G = (int)iGrid, hg = (G + 1)/2;
    const size_t half = (size_t)iHalfGridNum;
    const bool salt = fZero < std::abs(fIonStrength);
    if (phimap1.size() < half || phimap2.size() < half) return false;
    const int nDiel = (int)iDielecBndyOdd, nDielEven = (int)iDielecBndyEven;
    const int nChg = (int)iCrgedGridSum, nChgEven = (int)iCrgedGridEven;

    // two charges on one grid point would race in add_charge: leave such a case to the CPU
    {
        std::vector<char> seen(half + 2, 0);
        for (int n = 0; n < nChg; n++)
        {
            if (n == nChgEven) std::fill(seen.begin(), seen.end(), 0);
            const long p = (long)prgiCrgPose[n];
            if (p < 1 || p > (long)half) return false;
            if (seen[p]) return false;
            seen[p] = 1;
        }
    }

    const auto t0 = std::chrono::steady_clock::now();

    // ---- host copies in the layout of the kernels
    std::vector<int> hSta1(sta1.begin(), sta1.end()), hFi1(fi1.begin(), fi1.end()), hSta2(sta2.begin(), sta2.end()), hFi2(fi2.begin(), fi2.end());
    std::vector<int> hIdpos((size_t)nDiel), hCpos((size_t)nChg);
    for (int n = 0; n < nDiel; n++) hIdpos[n] = (int)prgiBndyDielecIndex[n];
    for (int n = 0; n < nChg; n++) hCpos[n] = (int)prgiCrgPose[n];
    std::vector<delphi_real> hDb((size_t)nDiel*6), hCval(prgfCrgValA.begin(), prgfCrgValA.begin() + nChg);
    for (int n = 0; n < nDiel; n++) for (int k = 0; k < 6; k++) hDb[(size_t)n*6 + k] = prgfBndyDielec[n][k];
    size_t maxLen = 1;
    for (int n = 1; n < G - 1; n++)
    {
        maxLen = std::max(maxLen, (size_t)(hFi1[n] - hSta1[n] + 1));
        maxLen = std::max(maxLen, (size_t)(hFi2[n] - hSta2[n] + 1));
    }
    const int cntOdd  = (G*(G-1) - 2)/2 - (G + 1)/2 + 1;
    const int cntEven = (G*(G-1) - 1)/2 - (G + 2)/2 + 1;
    if (cntOdd < 0 || cntEven < 0 || (int)bndx1.size() < cntOdd || (int)bndx2.size() < cntOdd || (int)bndx3.size() < cntEven || (int)bndx4.size() < cntEven) return false;
    const int iy0Odd = G*(G + 1)/2 - G + 1, iy0Even = G*(G + 2)/2 - G + 1;

    // ---- device buffers
    std::vector<cl_mem> bufs;
    bool failed = false;
    auto mk = [&](size_t bytes, const void* host, bool ro) -> cl_mem { cl_mem m = R.buffer(bytes, host, ro); if (!m) failed = true; else bufs.push_back(m); return m; };
    const size_t nGroups = 64;
    std::vector<delphi_real> hPmax(nGroups), hPsum(nGroups);
    cl_mem dPhi1 = mk(half*sizeof(delphi_real), phimap1.data(), false), dPhi2 = mk(half*sizeof(delphi_real), phimap2.data(), false);
    cl_mem dSnap = mk(half*sizeof(delphi_real), phimap2.data(), false);
    cl_mem dSf1 = mk(salt ? half*sizeof(delphi_real) : 0, salt ? prgfSaltMap1.data() : 0, true), dSf2 = mk(salt ? half*sizeof(delphi_real) : 0, salt ? prgfSaltMap2.data() : 0, true);
    cl_mem dSta1 = mk(hSta1.size()*4, hSta1.data(), true), dFi1 = mk(hFi1.size()*4, hFi1.data(), true);
    cl_mem dSta2 = mk(hSta2.size()*4, hSta2.data(), true), dFi2 = mk(hFi2.size()*4, hFi2.data(), true);
    cl_mem dIdpos = mk(hIdpos.size()*4, hIdpos.data(), true), dDb = mk(hDb.size()*sizeof(delphi_real), hDb.data(), true);
    cl_mem dBx1 = mk(cntOdd*sizeof(delphi_real), bndx1.data(), true), dBx2 = mk(cntOdd*sizeof(delphi_real), bndx2.data(), true);
    cl_mem dBx3 = mk(cntEven*sizeof(delphi_real), bndx3.data(), true), dBx4 = mk(cntEven*sizeof(delphi_real), bndx4.data(), true);
    cl_mem dCpos = mk(hCpos.size()*4, hCpos.data(), true), dCval = mk(hCval.size()*sizeof(delphi_real), hCval.data(), true);
    cl_mem dPmax = mk(nGroups*sizeof(delphi_real), hPmax.data(), false), dPsum = mk(nGroups*sizeof(delphi_real), hPsum.data(), false);
    auto release = [&]() { for (size_t i = 0; i < bufs.size(); i++) A.ReleaseMemObject(bufs[i]); bufs.clear(); };
    if (failed) { release(); return false; }

    // ---- fixed kernel arguments
    cl_int err = 0;
    auto arg = [&](cl_kernel k, cl_uint i, size_t sz, const void* v) { err |= A.SetKernelArg(k, i, sz, v); };
    const delphi_real vOm1 = om1, vSixth = sixth;
    const int vSalt = salt ? 1 : 0, vLat1 = (int)lat1, vLat2 = (int)lat2, vLong1 = (int)long1, vLong2 = (int)long2;
    const int n0Zero = 0, nHalfGrid = (int)half;
    for (int odd = 0; odd < 2; odd++)
    {
        cl_kernel k = odd ? R.kSweepOdd : R.kSweepEven;
        cl_mem a = odd ? dPhi1 : dPhi2, b = odd ? dPhi2 : dPhi1, sf = odd ? dSf1 : dSf2, st = odd ? dSta1 : dSta2, fn = odd ? dFi1 : dFi2;
        arg(k, 0, sizeof a, &a); arg(k, 1, sizeof b, &b); arg(k, 2, sizeof sf, &sf); arg(k, 3, sizeof st, &st); arg(k, 4, sizeof fn, &fn);
        arg(k, 5, sizeof vOm1, &vOm1); arg(k, 6, sizeof vSixth, &vSixth); arg(k, 7, 4, &vSalt);
        arg(k, 8, 4, &vLat1); arg(k, 9, 4, &vLat2); arg(k, 10, 4, &vLong1); arg(k, 11, 4, &vLong2);

        cl_kernel kc = odd ? R.kCorrOdd : R.kCorrEven;
        const int n0 = odd ? 0 : nDielEven;
        arg(kc, 0, sizeof a, &a); arg(kc, 1, sizeof b, &b); arg(kc, 2, sizeof dIdpos, &dIdpos); arg(kc, 3, sizeof dDb, &dDb);
        arg(kc, 4, 4, &n0); arg(kc, 5, 4, &vLat1); arg(kc, 6, 4, &vLat2); arg(kc, 7, 4, &vLong1); arg(kc, 8, 4, &vLong2);
    }
    // reset_faces and add_charge are shared by both halves, so their arguments are set before every launch (cheap)
    auto launch2D = [&](cl_kernel k, size_t g0, size_t g1) { size_t g[2] = { g0, g1 }; err |= A.EnqueueNDRangeKernel(R.q, k, 2, 0, g, 0, 0, 0, 0); };
    auto launch1D = [&](cl_kernel k, size_t n) { if (n) { size_t g = n; err |= A.EnqueueNDRangeKernel(R.q, k, 1, 0, &g, 0, 0, 0, 0); } };
    const size_t g0 = (maxLen + 63)/64*64, g1 = (size_t)(G - 2);
    auto halfSweep = [&](bool odd)
    {
        launch2D(odd ? R.kSweepOdd : R.kSweepEven, g0, g1);
        launch1D(odd ? R.kCorrOdd : R.kCorrEven, (size_t)(odd ? nDielEven : nDiel - nDielEven));
        cl_mem a = odd ? dPhi1 : dPhi2, bx1 = odd ? dBx1 : dBx3, bx2 = odd ? dBx2 : dBx4;
        const int iy0 = odd ? iy0Odd : iy0Even;
        arg(R.kReset, 0, sizeof a, &a); arg(R.kReset, 1, sizeof bx1, &bx1); arg(R.kReset, 2, sizeof bx2, &bx2); arg(R.kReset, 3, 4, &iy0); arg(R.kReset, 4, 4, &G); arg(R.kReset, 5, 4, &hg);
        launch1D(R.kReset, (size_t)(odd ? cntOdd : cntEven));
        const int c0 = odd ? 0 : nChgEven;
        arg(R.kCharge, 0, sizeof a, &a); arg(R.kCharge, 1, sizeof dCpos, &dCpos); arg(R.kCharge, 2, sizeof dCval, &dCval); arg(R.kCharge, 3, 4, &c0);
        launch1D(R.kCharge, (size_t)(odd ? nChgEven : nChg - nChgEven));
    };
    arg(R.kReduce, 0, sizeof dSnap, &dSnap); arg(R.kReduce, 1, sizeof dPhi2, &dPhi2); arg(R.kReduce, 2, 4, &nHalfGrid);
    arg(R.kReduce, 3, sizeof dPmax, &dPmax); arg(R.kReduce, 4, sizeof dPsum, &dPsum);
    arg(R.kReduce, 5, R.workGroup*sizeof(delphi_real), 0); arg(R.kReduce, 6, R.workGroup*sizeof(delphi_real), 0);
    if (err) { release(); return false; }
    (void)n0Zero;

    // ---- the loop of itit(), with the same bookkeeping
    const int interval = iIterateInterval > 0 ? iIterateInterval : 1;
    int itr = 1, ires = 0;
    do
    {
        halfSweep(true);
        halfSweep(false);
        if (interval - 1 == itr%interval) err |= A.EnqueueCopyBuffer(R.q, dPhi2, dSnap, 0, 0, half*sizeof(delphi_real), 0, 0, 0);
        if (0 == itr%interval)
        {
            size_t gl = nGroups*R.workGroup, ll = R.workGroup;
            err |= A.EnqueueNDRangeKernel(R.q, R.kReduce, 1, 0, &gl, &ll, 0, 0, 0);
            err |= A.EnqueueReadBuffer(R.q, dPmax, 1, 0, nGroups*sizeof(delphi_real), hPmax.data(), 0, 0, 0);
            err |= A.EnqueueReadBuffer(R.q, dPsum, 1, 0, nGroups*sizeof(delphi_real), hPsum.data(), 0, 0, 0);
            if (err) { release(); return false; }
            double rn = 0.0, mx = 0.0;
            for (size_t g = 0; g < nGroups; g++) { rn += (double)hPsum[g]; if ((double)hPmax[g] > mx) mx = (double)hPmax[g]; }
            const delphi_real rmsch = (delphi_real)std::sqrt(rn/((double)(G - 2)*(G - 2)*(G - 2)));
            const delphi_real rmxch = (delphi_real)mx;
            if (fRmsc > rmsch || fMaxc > rmxch) ires = 1;
        }
        if (err) { release(); return false; }

        itr++;
        if (itr > delphi::DELPHI_MAX_ITERATIONS)
        {
            cerr << " WARNING: SOR stopped after " << delphi::DELPHI_MAX_ITERATIONS << " iterations without reaching the convergence threshold." << endl;
            break;
        }
        if (1.0e-7 > maxres) { if (iLinIterateNum >= itr && 0 == ires) continue; else break; }
        else                 { if ((iLinIterateNum >= itr || bAutoConverge) && 0 == ires) continue; else break; }
    }
    while (true);

    // ---- results back to the host: into temporaries first, so a failed read cannot leave the host arrays half written
    std::vector<delphi_real> out1(half), out2(half);
    err |= A.EnqueueReadBuffer(R.q, dPhi1, 1, 0, half*sizeof(delphi_real), out1.data(), 0, 0, 0);
    err |= A.EnqueueReadBuffer(R.q, dPhi2, 1, 0, half*sizeof(delphi_real), out2.data(), 0, 0, 0);
    A.Finish(R.q);
    release();
    if (err) return false;
    std::copy(out1.begin(), out1.end(), phimap1.begin());
    std::copy(out2.begin(), out2.end(), phimap2.begin());
    itrOut = itr;

    if (std::getenv("DELPHI_GPU_VERBOSE"))
        std::fprintf(stderr, " GPU solve: grid %d, %d iterations, %.1f ms\n", G, itr, 1e3*std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return true;
}
