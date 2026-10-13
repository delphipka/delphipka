// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: new file.

//
// Minimal OpenCL runtime for the optional GPU path of the SOR solver (solver_fastSOR_gpu.cpp). Linux, macOS and Windows.
//
// OpenCL is vendor-neutral (NVIDIA, AMD, Intel, and Apple through Mesa's Rusticl). It is loaded at run time with dlopen, so
// building needs no OpenCL headers or libraries, and a machine without a usable GPU simply keeps running on the CPU.
//   DELPHI_GPU=0            never use the GPU
//   DELPHI_GPU=1            use it even for small grids (testing)
//   DELPHI_GPU_MIN_GRID=N   smallest grid edge sent to the GPU (default GPU_MIN_GRID below)
//   DELPHI_GPU_DEVICE=k     use the k-th GPU instead of the fastest one
//   DELPHI_GPU_VERBOSE=1    report every solve done on the GPU
#ifndef DELPHI_GPU_CL_H_
#define DELPHI_GPU_CL_H_

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// Smallest grid edge for which a GPU solve is worth its launch overhead (measured on an Apple M2: the CPU wins below ~100).
const int GPU_MIN_GRID = 101;

namespace delphicl
{
// dlopen/dlsym on Linux and macOS, LoadLibrary/GetProcAddress on Windows
inline void* open_library(const char* name)
{
#if defined(_WIN32)
    return (void*)LoadLibraryA(name);
#else
    return dlopen(name, RTLD_NOW | RTLD_LOCAL);
#endif
}
inline void* find_symbol(void* lib, const char* name)
{
#if defined(_WIN32)
    return (void*)GetProcAddress((HMODULE)lib, name);
#else
    return dlsym(lib, name);
#endif
}

typedef int32_t cl_int;  typedef uint32_t cl_uint;  typedef uint64_t cl_ulong;
typedef struct _p*  cl_platform_id;  typedef struct _d*  cl_device_id;  typedef struct _c*  cl_context;
typedef struct _q*  cl_command_queue;  typedef struct _m*  cl_mem;  typedef struct _pr* cl_program;
typedef struct _k*  cl_kernel;  typedef struct _e*  cl_event;

enum { DEVICE_TYPE_GPU = 4, MEM_READ_WRITE = 1, MEM_READ_ONLY = 4, MEM_COPY_HOST_PTR = 32,
       DEVICE_NAME = 0x102B, DEVICE_MAX_COMPUTE_UNITS = 0x1002, DEVICE_MAX_CLOCK = 0x100C, DEVICE_DOUBLE_FP = 0x1032,
       DEVICE_MAX_WORK_GROUP = 0x1004, PROGRAM_BUILD_LOG = 0x1183 };

struct Api
{
    cl_int (*GetPlatformIDs)(cl_uint, cl_platform_id*, cl_uint*);
    cl_int (*GetDeviceIDs)(cl_platform_id, cl_ulong, cl_uint, cl_device_id*, cl_uint*);
    cl_int (*GetDeviceInfo)(cl_device_id, cl_uint, size_t, void*, size_t*);
    cl_context (*CreateContext)(const intptr_t*, cl_uint, const cl_device_id*, void (*)(const char*, const void*, size_t, void*), void*, cl_int*);
    cl_command_queue (*CreateCommandQueueWithProperties)(cl_context, cl_device_id, const cl_ulong*, cl_int*);
    cl_command_queue (*CreateCommandQueue)(cl_context, cl_device_id, cl_ulong, cl_int*);
    cl_mem (*CreateBuffer)(cl_context, cl_ulong, size_t, void*, cl_int*);
    cl_int (*ReleaseMemObject)(cl_mem);
    cl_program (*CreateProgramWithSource)(cl_context, cl_uint, const char**, const size_t*, cl_int*);
    cl_int (*BuildProgram)(cl_program, cl_uint, const cl_device_id*, const char*, void (*)(cl_program, void*), void*);
    cl_int (*GetProgramBuildInfo)(cl_program, cl_device_id, cl_uint, size_t, void*, size_t*);
    cl_kernel (*CreateKernel)(cl_program, const char*, cl_int*);
    cl_int (*ReleaseKernel)(cl_kernel);
    cl_int (*SetKernelArg)(cl_kernel, cl_uint, size_t, const void*);
    cl_int (*EnqueueNDRangeKernel)(cl_command_queue, cl_kernel, cl_uint, const size_t*, const size_t*, const size_t*, cl_uint, const cl_event*, cl_event*);
    cl_int (*EnqueueReadBuffer)(cl_command_queue, cl_mem, cl_uint, size_t, size_t, void*, cl_uint, const cl_event*, cl_event*);
    cl_int (*EnqueueCopyBuffer)(cl_command_queue, cl_mem, cl_mem, size_t, size_t, size_t, cl_uint, const cl_event*, cl_event*);
    cl_int (*Finish)(cl_command_queue);
};

// OpenCL C source of the kernels. REAL is float (default build) or double; see solver_fastSOR_gpu.cpp for what each one reproduces.
static const char* KERNEL_SOURCE = R"CLC(
#ifdef USE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

// main stencil of one half sweep over the planes n = 1..G-2, point ix in [sta[n], fin[n]]  (odd half: a = phimap1, b = phimap2)
__kernel void sweep_odd(__global REAL* a, __global const REAL* b, __global const REAL* sf, __global const int* sta, __global const int* fin,
                        const REAL om1, const REAL sixth, const int salt, const int lat1, const int lat2, const int long1, const int long2)
{
    int n = 1 + (int)get_global_id(1);
    int ix = sta[n] + (int)get_global_id(0);
    if (ix > fin[n]) return;
    int j = ix - 1;
    REAL t1 = b[j] + b[j-1];
    REAL t2 = b[j+lat1] + b[j-lat2];
    REAL t3 = b[j+long1] + b[j-long2];
    REAL s = salt ? sf[j] : sixth;
    a[j] = a[j]*om1 + (t1+t2+t3)*s;
}

// even half: a = phimap2, b = phimap1
__kernel void sweep_even(__global REAL* a, __global const REAL* b, __global const REAL* sf, __global const int* sta, __global const int* fin,
                         const REAL om1, const REAL sixth, const int salt, const int lat1, const int lat2, const int long1, const int long2)
{
    int n = 1 + (int)get_global_id(1);
    int ix = sta[n] + (int)get_global_id(0);
    if (ix > fin[n]) return;
    int j = ix - 1;
    REAL t1 = b[j] + b[j+1];
    REAL t2 = b[j+lat2] + b[j-lat1];
    REAL t3 = b[j+long2] + b[j-long1];
    REAL s = salt ? sf[j] : sixth;
    a[j] = a[j]*om1 + (t1+t2+t3)*s;
}

// dielectric boundary points, entries n0 .. n0+count-1 of the list (odd half: the first icount2a, even half: the rest)
__kernel void corr_odd(__global REAL* a, __global const REAL* b, __global const int* idpos, __global const REAL* db,
                       const int n0, const int lat1, const int lat2, const int long1, const int long2)
{
    int n = n0 + (int)get_global_id(0);
    int ix = idpos[n];
    __global const REAL* d = db + 6*n;
    REAL t1 = b[ix-2]*d[0]       + b[ix-1]*d[1];
    REAL t2 = b[ix-1-lat2]*d[2]  + b[ix-1+lat1]*d[3];
    REAL t3 = b[ix-1-long2]*d[4] + b[ix-1+long1]*d[5];
    a[ix-1] += t1 + t2 + t3;
}

__kernel void corr_even(__global REAL* a, __global const REAL* b, __global const int* idpos, __global const REAL* db,
                        const int n0, const int lat1, const int lat2, const int long1, const int long2)
{
    int n = n0 + (int)get_global_id(0);
    int ix = idpos[n];
    __global const REAL* d = db + 6*n;
    REAL t1 = b[ix-1]*d[0]       + b[ix]*d[1];
    REAL t2 = b[ix-1-lat1]*d[2]  + b[ix-1+lat2]*d[3];
    REAL t3 = b[ix-1-long1]*d[4] + b[ix-1+long2]*d[5];
    a[ix-1] += t1 + t2 + t3;
}

// restore the box-face values the plane loops overwrote:  iy = iy0 + (n+1)*G
__kernel void reset_faces(__global REAL* a, __global const REAL* bx1, __global const REAL* bx2, const int iy0, const int G, const int hg)
{
    int n = (int)get_global_id(0);
    int iy = iy0 + (n+1)*G;
    a[iy-1] = bx1[n];
    a[iy + (hg-1) - 1] = bx2[n];
}

// charge terms, entries n0 .. n0+count-1 (grid points are distinct, checked on the host)
__kernel void add_charge(__global REAL* a, __global const int* pos, __global const REAL* val, const int n0)
{
    int n = n0 + (int)get_global_id(0);
    a[pos[n]-1] += val[n];
}

// max and sum of squares of (snap - cur) over i = 1 .. n-1, one partial result per work-group
__kernel void diff_reduce(__global const REAL* snap, __global const REAL* cur, const int n, __global REAL* pmax, __global REAL* psum,
                          __local REAL* lmax, __local REAL* lsum)
{
    int gid = (int)get_global_id(0), gs = (int)get_global_size(0), lid = (int)get_local_id(0);
    REAL mx = 0, sm = 0;
    for (int i = 1 + gid; i < n; i += gs) {
        REAL d = snap[i] - cur[i];
        REAL ad = d < 0 ? -d : d;
        mx = ad > mx ? ad : mx;
        sm += d*d;
    }
    lmax[lid] = mx; lsum[lid] = sm;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = (int)get_local_size(0)/2; s > 0; s >>= 1) {
        if (lid < s) { lmax[lid] = lmax[lid] > lmax[lid+s] ? lmax[lid] : lmax[lid+s]; lsum[lid] += lsum[lid+s]; }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (lid == 0) { pmax[get_group_id(0)] = lmax[0]; psum[get_group_id(0)] = lsum[0]; }
}
)CLC";

// One process-wide context on the chosen GPU, created on first use. ok() is false when OpenCL or a suitable GPU is missing.
class Runtime
{
public:
    // fp64: true when the solver's reals are double; the first call decides (the precision is fixed per build)
    static Runtime& get(bool fp64) { static Runtime r(fp64); return r; }
    bool ok() const { return ok_; }
    Api api;
    cl_context ctx; cl_command_queue q; cl_device_id dev;
    std::string deviceName;
    size_t workGroup;                 // work-group size used by the reduction (power of two)
    cl_kernel kSweepOdd, kSweepEven, kCorrOdd, kCorrEven, kReset, kCharge, kReduce;

    // makes a read-only / read-write device buffer initialised from the host
    cl_mem buffer(size_t bytes, const void* host, bool readOnly) const
    {
        cl_int err = 0;
        cl_mem m = api.CreateBuffer(ctx, (readOnly ? MEM_READ_ONLY : MEM_READ_WRITE) | MEM_COPY_HOST_PTR, bytes ? bytes : 4, const_cast<void*>(host ? host : (const void*)&dummy_), &err);
        return err ? 0 : m;
    }

private:
    bool ok_;
    int dummy_;
    void* lib_;
    explicit Runtime(bool fp64) : ok_(false), dummy_(0), lib_(0) { ok_ = init(fp64); }

    template <class F> bool sym(F& f, const char* name, bool required = true)
    {
        f = reinterpret_cast<F>(find_symbol(lib_, name));
        return f != 0 || !required;
    }

    bool init(bool fp64)
    {
        const char* off = std::getenv("DELPHI_GPU");
        if (off && '0' == off[0]) return false;
#if defined(_WIN32)
        const char* names[] = { "OpenCL.dll", 0 };
#elif defined(__APPLE__)
        const char* names[] = { "/System/Library/Frameworks/OpenCL.framework/OpenCL", "libOpenCL.dylib", 0 };
#else
        const char* names[] = { "libOpenCL.so.1", "libOpenCL.so", 0 };
#endif
        for (int i = 0; names[i] && !lib_; i++) lib_ = open_library(names[i]);
        if (!lib_) return false;
        if (!(sym(api.GetPlatformIDs, "clGetPlatformIDs") && sym(api.GetDeviceIDs, "clGetDeviceIDs") && sym(api.GetDeviceInfo, "clGetDeviceInfo")
              && sym(api.CreateContext, "clCreateContext") && sym(api.CreateBuffer, "clCreateBuffer") && sym(api.ReleaseMemObject, "clReleaseMemObject")
              && sym(api.CreateProgramWithSource, "clCreateProgramWithSource") && sym(api.BuildProgram, "clBuildProgram")
              && sym(api.GetProgramBuildInfo, "clGetProgramBuildInfo") && sym(api.CreateKernel, "clCreateKernel") && sym(api.ReleaseKernel, "clReleaseKernel")
              && sym(api.SetKernelArg, "clSetKernelArg") && sym(api.EnqueueNDRangeKernel, "clEnqueueNDRangeKernel")
              && sym(api.EnqueueReadBuffer, "clEnqueueReadBuffer") && sym(api.EnqueueCopyBuffer, "clEnqueueCopyBuffer") && sym(api.Finish, "clFinish")))
            return false;
        sym(api.CreateCommandQueueWithProperties, "clCreateCommandQueueWithProperties", false);
        sym(api.CreateCommandQueue, "clCreateCommandQueue", false);

        // the fastest GPU (compute units x clock), or the one named by DELPHI_GPU_DEVICE
        std::vector<cl_device_id> gpus;
        cl_platform_id plats[16]; cl_uint np = 0;
        if (api.GetPlatformIDs(16, plats, &np) || 0 == np) return false;
        for (cl_uint p = 0; p < np; p++)
        {
            cl_device_id ds[16]; cl_uint nd = 0;
            if (api.GetDeviceIDs(plats[p], DEVICE_TYPE_GPU, 16, ds, &nd) == 0)
                for (cl_uint d = 0; d < nd; d++) gpus.push_back(ds[d]);
        }
        if (gpus.empty()) return false;
        size_t pick = 0; double best = -1;
        const char* want = std::getenv("DELPHI_GPU_DEVICE");
        if (want) pick = (size_t)std::atoi(want) < gpus.size() ? (size_t)std::atoi(want) : 0;
        else for (size_t i = 0; i < gpus.size(); i++)
        {
            cl_uint cu = 0, mhz = 0; api.GetDeviceInfo(gpus[i], DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, 0); api.GetDeviceInfo(gpus[i], DEVICE_MAX_CLOCK, sizeof mhz, &mhz, 0);
            if ((double)cu*(mhz ? mhz : 1) > best) { best = (double)cu*(mhz ? mhz : 1); pick = i; }
        }
        dev = gpus[pick];
        char nm[256] = {0}; api.GetDeviceInfo(dev, DEVICE_NAME, sizeof nm, nm, 0); deviceName = nm;

        if (fp64) { cl_ulong cfg = 0; api.GetDeviceInfo(dev, DEVICE_DOUBLE_FP, sizeof cfg, &cfg, 0); if (!cfg) return false; }

        cl_int err = 0;
        ctx = api.CreateContext(0, 1, &dev, 0, 0, &err); if (err || !ctx) return false;
        q = 0;
        if (api.CreateCommandQueueWithProperties) q = api.CreateCommandQueueWithProperties(ctx, dev, 0, &err);
        if ((err || !q) && api.CreateCommandQueue) { err = 0; q = api.CreateCommandQueue(ctx, dev, 0, &err); }
        if (err || !q) return false;

        std::string src = std::string("#define REAL ") + (fp64 ? "double\n#define USE_FP64 1\n" : "float\n") + KERNEL_SOURCE;
        const char* sp = src.c_str(); size_t sl = src.size();
        cl_program prog = api.CreateProgramWithSource(ctx, 1, &sp, &sl, &err); if (err || !prog) return false;
        if (api.BuildProgram(prog, 1, &dev, "", 0, 0))
        {
            char log[2048] = {0}; api.GetProgramBuildInfo(prog, dev, PROGRAM_BUILD_LOG, sizeof log - 1, log, 0);
            std::fprintf(stderr, " GPU: OpenCL build failed on %s, staying on the CPU:\n%s\n", deviceName.c_str(), log);
            return false;
        }
        kSweepOdd = api.CreateKernel(prog, "sweep_odd", &err);   if (err) return false;
        kSweepEven = api.CreateKernel(prog, "sweep_even", &err); if (err) return false;
        kCorrOdd = api.CreateKernel(prog, "corr_odd", &err);     if (err) return false;
        kCorrEven = api.CreateKernel(prog, "corr_even", &err);   if (err) return false;
        kReset = api.CreateKernel(prog, "reset_faces", &err);    if (err) return false;
        kCharge = api.CreateKernel(prog, "add_charge", &err);    if (err) return false;
        kReduce = api.CreateKernel(prog, "diff_reduce", &err);   if (err) return false;
        size_t maxwg = 64; api.GetDeviceInfo(dev, DEVICE_MAX_WORK_GROUP, sizeof maxwg, &maxwg, 0);
        workGroup = 1; while (workGroup*2 <= maxwg && workGroup < 128) workGroup *= 2;
        if (std::getenv("DELPHI_GPU_VERBOSE")) std::fprintf(stderr, " GPU: using %s (OpenCL, %s)\n", deviceName.c_str(), fp64 ? "double" : "single precision");
        return true;
    }
};
} // namespace delphicl

#endif
