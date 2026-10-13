// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: new file.

//
// Pins each MPI rank to its own cores from inside the program, so no shell wrapper (taskset, numactl) is needed and the same
// command works on Linux and Windows:  mpiexec -n 4 delphiPKa run.prm   with the environment below.
//   DELPHI_FIRST_CORE=f   rank r gets the cores f + r*T .. f + r*T + T-1   (unset: no pinning, the OS/MPI decide)
//   DELPHI_THREADS=T      cores per rank (default 1); also the number of OpenMP threads the launcher asks for
//   DELPHI_PIN=0          never pin
//   DELPHI_PIN_VERBOSE=1  report the cores each rank got
// Linux: sched_setaffinity. Windows: SetProcessAffinityMask (the first 64 logical processors). macOS has no hard thread affinity,
// so the request is ignored there. On machines with fast and slow cores give FIRST_CORE the first fast core.
#ifndef DELPHI_CPU_PIN_H_
#define DELPHI_CPU_PIN_H_

#include <cstdio>
#include <cstdlib>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <sched.h>
#endif

inline void delphi_pin_rank(int rank)
{
    const char* first = std::getenv("DELPHI_FIRST_CORE");
    const char* off   = std::getenv("DELPHI_PIN");
    if (!first || (off && '0' == off[0])) return;
    const char* th = std::getenv("DELPHI_THREADS");
    int t = th ? std::atoi(th) : 1;
    if (t < 1) t = 1;
    const int lo = std::atoi(first) + rank*t, hi = lo + t - 1;
    int rc = 0;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    for (int c = lo; c <= hi && c < CPU_SETSIZE; c++) CPU_SET(c, &set);
    rc = sched_setaffinity(0, sizeof set, &set);
#elif defined(_WIN32)
    DWORD_PTR mask = 0;
    for (int c = lo; c <= hi && c < 64; c++) mask |= ((DWORD_PTR)1) << c;
    rc = (mask && SetProcessAffinityMask(GetCurrentProcess(), mask)) ? 0 : -1;
#else
    rc = 0;   // no hard affinity on this platform
#endif
    if (rc != 0)
        std::fprintf(stderr, " WARNING: rank %d could not be pinned to cores %d-%d (does the machine have them?), continuing unpinned.\n", rank, lo, hi);
    else if (std::getenv("DELPHI_PIN_VERBOSE"))
        std::fprintf(stderr, " rank %d pinned to cores %d-%d\n", rank, lo, hi);
}

#endif
