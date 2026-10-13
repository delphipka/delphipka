// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: new file. Uninitialised buffers; gepsmp only when the second Gaussian run needs it.

//
// Reuse of the products of CDelphiSpace::setGaussian(). The charged and the neutral state of a titratable site are
// solved on the same grids around the same atoms and differ only in a few charges, which setGaussian() never reads,
// so its dielectric maps and boundary-point list are identical and the second state can take them from here.
// An entry is only used when every input of setGaussian() is bitwise equal (the key below), so a hit changes nothing.
//   DELPHI_GAUSS_CACHE=0        switches the cache off
//   DELPHI_GAUSS_CACHE_CHECK=1  on a hit also recomputes and aborts if the restored state differs (test aid)
#ifndef DELPHI_GAUSS_CACHE_H_
#define DELPHI_GAUSS_CACHE_H_

#include <cstddef>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../interface/environment.h"
#include "../misc/misc_grid.h"

const size_t GAUSS_CACHE_ENTRIES     = 4;                    // one site's four solves (base, two focusings, residue alone)
const size_t GAUSS_CACHE_ENTRY_BYTES = 96ul*1024*1024;       // larger states are not kept: for a 175^3 grid (265 MB) copying
                                                             // them in and out costs more than setGaussian() and, with several
                                                             // ranks, the memory makes the whole run slower (3HDD: 4 ranks
                                                             // 252 s with the cache, 135 s without)
const size_t GAUSS_CACHE_BYTES       = 256ul*1024*1024;      // upper bound on what is kept per process

// Uninitialised array: the maps are overwritten in full on capture, so the zero fill of a std::vector (about 1.3 s on a 62-site run) is waste.
template <class T> struct RawBuf
{
    std::unique_ptr<T[]> p;
    size_t n;
    RawBuf() : n(0) {}
    void alloc(size_t k) { p.reset(new T[k]); n = k; }
    T* data() { return p.get(); }
    const T* data() const { return p.get(); }
    size_t size() const { return n; }
};

// What setGaussian() leaves behind that later code reads: gepsmp2, idebmap, the boundary-point list and its length. gepsmp (the
// density) is only read back by the second Gaussian run (logs), so it is kept only then (hasG1).
struct GaussState
{
    RawBuf< SGrid<delphi_real> > g1, g2;                // gepsmp (only if hasG1), gepsmp2, over indices 0..iGrid
    RawBuf<char> deb;                                   // idebmap over indices 0..iGrid
    std::vector< SGrid<delphi_integer> > ibgrd;         // boundary points
    delphi_integer ibound;                              // iBoundNum
    bool hasG1;

    GaussState() : ibound(0), hasG1(false) {}
    size_t bytes() const
    {
        return (g1.size() + g2.size())*sizeof(SGrid<delphi_real>) + deb.size() + ibgrd.size()*sizeof(SGrid<delphi_integer>);
    }
};

class GaussCache
{
public:
    struct Entry { std::string key; GaussState state; Entry() {} Entry(Entry&& o) : key(std::move(o.key)), state(std::move(o.state)) {} Entry& operator=(Entry&& o) { key = std::move(o.key); state = std::move(o.state); return *this; } };

    static GaussCache& instance() { static GaussCache c; return c; }

    const GaussState* find(const std::string& key) const
    {
        for (const Entry& e : entries) if (e.key == key) return &e.state;
        return NULL;
    }

    void store(const std::string& key, GaussState& state)
    {
        if (state.bytes() > GAUSS_CACHE_ENTRY_BYTES) return;
        entries.push_back(Entry());
        entries.back().key = key;
        entries.back().state = std::move(state);
        size_t total = 0;
        for (const Entry& e : entries) total += e.state.bytes();
        while (!entries.empty() && (entries.size() > GAUSS_CACHE_ENTRIES || total > GAUSS_CACHE_BYTES))
        {
            total -= entries.front().state.bytes();
            entries.pop_front();
        }
    }

    template <class T> static void add(std::string& key, const T& v) { key.append(reinterpret_cast<const char*>(&v), sizeof(T)); }

private:
    std::deque<Entry> entries;
};

#endif
