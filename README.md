 

DelphiPKA is a DelPhi-based C++ program, allowing to predict pKa's for ionizable groups in proteins, RNA and DNA.
Unique approach stems from:
    a) Use gaussian-based smooth function to mimic conformational changes associated with ionization changes.
    b) Calculate the electrostatic energy without defining molecular surface.

For questions and help, visit http://compbio.clemson.edu/forum/ or email to delphi@g.clemson.edu

User Manual can be found [here ](http://compbio.clemson.edu/pka_webserver/assets/manual/DelPhiPKa_User_Manual.html)

If the use of DelPhiPKa results to a publication, please cite the following references:

*1.   Lin Wang, Lin Li, and Emil Alexov. "pKa predictions for proteins, RNAs and DNAs with the Gaussian dielectric function using DelPhiPKa." Proteins. (2015) Sep 26, doi: 10.1002/prot.24935*

*2.   Lin Wang, Min Zhang, and Emil Alexov. "DelPhiPKa Web Server: Predicting pKa of proteins, RNAs and DNAs." Bioinformatics. (2015) Oct 29, doi: 10.1093/bioinformatics/btv607*

## Building and running (branch nucleic-acids)

Needs a C++11 compiler with OpenMP, Open MPI (`mpic++`), GSL and the Boost headers (`boost/lexical_cast.hpp`, `boost/any.hpp`).

On Fedora:

```bash
sudo dnf install gcc-c++ openmpi-devel gsl-devel boost-devel
export PATH=/usr/lib64/openmpi/bin:$PATH
make -j4
```

Run it with one process, or with several MPI ranks that each take a share of the titratable sites. `scripts/run_parallel.py` does
the second on Linux, macOS and Windows (Python 3, standard library only; it finds `mpiexec` or `mpirun` and recognises Open MPI,
MPICH, Intel MPI and Microsoft MPI):

```bash
./delphiPKa run.prm                                          # one process, one thread
python scripts/run_parallel.py -n 4 run.prm                  # 4 MPI ranks
python scripts/run_parallel.py -n 4 -t 2 run.prm             # 4 ranks, 2 OpenMP threads inside each solve
python scripts/run_parallel.py -n 4 -f 4 --dry-run run.prm   # only print the command
```

The options are `-n` ranks, `-t` threads per rank, `-f` first core, `-x` the binary (default `./delphiPKa`, `delphiPKa.exe` on
Windows). Each rank pins itself to its own cores from inside the program (`src/delphiPKa/cpu_pin.h`), so no `taskset`-style wrapper
is needed; the results are identical for every split. Without the script, the same thing is a plain MPI command with environment
variables:

```bash
DELPHI_FIRST_CORE=0 DELPHI_THREADS=1 OMP_NUM_THREADS=1 mpiexec -n 4 ./delphiPKa run.prm
```

```powershell
$env:DELPHI_FIRST_CORE = "0"; $env:DELPHI_THREADS = "1"; $env:OMP_NUM_THREADS = "1"
mpiexec -n 4 .\delphiPKa.exe run.prm
```

On a machine with fast and slow cores the launcher starts on the fast ones by itself (Linux reports the core capacities; elsewhere give `-f`, or `DELPHI_FIRST_CORE` without the launcher), so keep ranks x threads within the
number of fast cores. `DELPHI_PIN=0` switches pinning off and `DELPHI_PIN_VERBOSE=1` reports the cores each rank got.

The SOR solver's relaxation factor is now computed from the grid (it used to be fixed at 0.8 and left the solver under-relaxed, so the pKa values stopped well short of convergence). Converged values differ from older output by up to about 0.4-0.6 pKa units on nucleic acids and about 0.1 on proteins.

### Platforms

| Platform | Pinning | GPU (OpenCL) | Status |
|---|---|---|---|
| Linux | `sched_setaffinity` | `libOpenCL.so.1` | developed and tested (Fedora, aarch64) |
| Windows | `SetProcessAffinityMask`, first 64 logical processors | `OpenCL.dll` | written for it, **not tested** |
| macOS | none (the OS has no hard thread affinity) | OpenCL framework | written for it, **not tested** |

The build needs an MPI compiler wrapper (`mpic++`), OpenMP, GSL and the Boost headers on every platform. On Windows use MSYS2
(MinGW-w64) with MS-MPI or MPICH, for example `pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gsl mingw-w64-x86_64-boost make`
plus the MPI package, and expect to adjust the `Makefile` include and library paths. The launcher, the pinning and the GPU loader are
platform-neutral; the original DelPhiPKa sources have not been checked for other POSIX-only calls.

### Precision and speed options

The code is built in **single precision** by default (`delphi_real` = `float`; the authors' `SP` option, repaired). It halves the memory traffic of the Poisson-Boltzmann solves and agrees with the double-precision build to the printed 0.01 pKa units on every structure tested (1BNA, 1SDR, 1RNA, 6PTI, 1CRN, 3HDD). The double-precision build is still available:

```bash
make EXTRA=-DMX
```

In single precision a per-iteration change below about 5e-5 cannot be reached, so `MAXC` and `MAXC0` are raised to at least 1e-4 (with a warning) and the SOR loop stops after 20,000 iterations with a warning. Two `run.prm` lines control the whole-molecule **base run**, which supplies the boundary values of the focusing runs and takes about half of the time of a large system (the focusing runs stay at scales 2 and 4). Paste them into `run.prm` exactly like this, each on
its own line (the values shown are the defaults):

```text
Base grid scale (SCALE0)                : 1.0
Convergence threshold, base run (MAXC0) : 0.01
```



`SCALE0` is the number of grid points per Angstrom of the base run. `MAXC0` is the convergence threshold of the base run only; left out, it is 10 x `MAXC`, which is 0.09 for the default `MAXC` = 0.009. `run.prm` only takes `#` comments at the start of a line, so do not append comments to these lines, and keep the space before the colon (the reader drops the character in front of it). The base run therefore converges to 10 x `MAXC` by default, and tightening `MAXC` tightens the base run with it. To get the base run as tight as the rest, set `MAXC0` equal to `MAXC`:

```text
Convergence threshold, base run (MAXC0) : 0.001
```



Calibration against fully converged references (single precision, one core). 3HDD is a 1,800-atom protein-DNA complex with 62 sites:


| SCALE0  | MAXC0                   | 3HDD time  | rms error | worst site |
| --------- | ------------------------- | ------------ | ----------- | ------------ |
| 1.0     | = MAXC (old behaviour)  | 95.0 s     | 0.003     | 0.01       |
| **1.0** | **10 x MAXC (default)** | **71.6 s** | **0.018** | **0.13**   |
| 0.9     | 0.003                   | 62.4 s     | 0.035     | 0.19       |
| 0.9     | 0.01                    | 55.5 s     | 0.040     | 0.20       |
| 0.8     | = MAXC                  | 52.6 s     | 0.037     | 0.17       |
| 0.8     | 0.01                    | 44.7 s     | 0.045     | 0.17       |

**`MAXC` default.** The default is `MAXC` = 0.009 (it was 0.001): it saves about a quarter of the SOR iterations and costs about twice the error. Against fully converged double-precision references, with the current code on 4 ranks: 3HDD rms 0.037 / worst site 0.17 at 0.009, against 0.010 / 0.13 at 0.001; 1BNA rms 0.0076 / worst 0.023 at 0.009 against 0.0038 / 0.010 at 0.001. The accuracy breaks down quickly beyond that (1BNA at `MAXC` = 0.012: worst site 0.15). Put `Convergence threshold (MAXC) : 0.001` into `run.prm` when the last 0.01 matters.

With `MAXC` = 0.001, the smaller structures stay within the printed 0.01 of the converged answer (1BNA 0.000, 1SDR 0.004, 1RNA 0.003, 6PTI 0.003, 1CRN 0.005 rms). The error is in pKa units against the converged answer; the few sites that carry it are DNA bases (3-5 of 62 sites above 0.05 for `SCALE0` 0.9 and 0.8). Proteins are insensitive (6PTI: rms 0.004 at `SCALE0` 0.7, 0.006 at 0.5); nucleic acids are not (rms at `SCALE0` 0.9: 1BNA 0.006, 1RNA 0.011, 1SDR 0.033; at 0.7: 0.010, 0.096, 0.051; below 0.7 a DNA/RNA run degrades quickly). For comparison, the original fixed relaxation factor left every pKa about 0.1-0.4 away from convergence. Set `MAXC0` equal to `MAXC` when the last 0.01 matters.

### Optional GPU (OpenCL)

The iteration loop of the large Poisson-Boltzmann solves (grids of 101^3 and up) runs on a GPU when one is available, through OpenCL: it is vendor-neutral (NVIDIA, AMD, Intel, and Apple through Mesa's Rusticl) and is loaded at run time, so building needs no OpenCL headers and a machine without a usable GPU simply stays on the CPU. Everything for one solve (the plane stencil, the dielectric-boundary correction, the face reset, the charge terms and the convergence reduction) stays on the device between iterations. It is single precision only, or double on devices that report `cl_khr_fp64` when built with `EXTRA=-DMX`.

Environment variables:

```bash
DELPHI_GPU=0            # never use the GPU
DELPHI_GPU=1            # use it for every grid, however small
DELPHI_GPU_MIN_GRID=N   # smallest grid edge sent to the GPU (default 101)
DELPHI_GPU_DEVICE=k     # use the k-th GPU instead of the fastest one
DELPHI_GPU_VERBOSE=1    # report each GPU solve
```

For example, to run on 4 ranks without the GPU:

```bash
DELPHI_GPU=0 python scripts/run_parallel.py -n 4 run.prm
```

With several MPI ranks the ranks share the GPU. Measured on an Apple M2 (the GPU shares the CPU's memory bandwidth, so the gain is
modest there; a discrete GPU should do much better), 3HDD with the default settings: 71.4 s -> 65.4 s on one rank, 39.9 s -> 34.3 s on
4 ranks, pKa values identical. Solves the kernels do not cover (periodic boundaries, the nonlinear solver, a grid-energy stop) and
any OpenCL error fall back to the CPU loop without changing the result.

### Code style

The C/C++ sources use 4-space indentation, no tabs and Allman braces; the rules are in `.clang-format`. `scripts/format.py` runs
`clang-format` on everything under `src/` (it works on Linux, macOS and Windows and needs `clang-format` on the `PATH`):

```bash
python scripts/format.py             # format the sources in place
python scripts/format.py --check     # only list the files that are not formatted (exit code 1 if there are any)
python scripts/format.py FILE ...    # only these files
```

The style changes whitespace (indentation, trailing spaces, `for(` to `for (`) and does not reflow lines, sort includes or touch
comments; the program is unchanged by it (a fully formatted build gives byte-identical results).
