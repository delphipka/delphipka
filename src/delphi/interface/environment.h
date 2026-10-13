// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: single precision is the default; limits for it.

/**
 * @file environment.h
 * @brief pre-processor marco's used for compiling the code
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

/**
 * preprocessor marco to choose data types for delphi_real and delphi_integer numbers \n
 * - MX -- double      &      int [mixed double/single precision, mostly used]
 * - SP -- float       &      int [single precision]
 * - DP -- double      & long int [double precision]
 * - LD -- long double & long int [extreme case,rarely used]
 */
#if !defined(SP) && !defined(MX) && !defined(DP) && !defined(LD)
#define SP   // default: single precision (build with  make EXTRA=-DMX  for double reals, as before)
#endif

/**
 * flag for additional outputs at run time
 */
//#define VERBOSE

/**
 * flag to control the precision of output numbers
 */
//#define DEVELOPER

/**
 * flag to embed delphicpp in mcce
 */
//#define MCCE

/**
 * flag to embed delphicpp in DelPhiPKa
 */
#define PRIME

/**
 * flag for OMP
 */
//#define PARALLEL_OMP

/**
 * Introduce namespace "delphi" here and used in getFunction and getStatement to avoid ambiguous reference
 * to data type "delphi_real" when compiling the code in Mac system
 */

namespace delphi
{
/**
 * Limits that keep a single-precision run finite: a float potential cannot get its per-iteration change below about
 * 5e-5 kT/e (1BNA; more for larger potentials), so MAXC is raised to at least DELPHI_MIN_MAXC, and the SOR loop stops
 * after DELPHI_MAX_ITERATIONS (normal runs need 50-1100) with a warning instead of iterating forever.
 */
#ifdef SP
const float DELPHI_MIN_MAXC = 1.0e-4f;
#else
const float DELPHI_MIN_MAXC = 0.0f;
#endif
const int DELPHI_MAX_ITERATIONS = 20000;

#ifdef SP
typedef int   delphi_integer;
typedef float delphi_real;
#endif

#ifdef MX
typedef int    delphi_integer;
typedef double delphi_real;
#endif

#ifdef DP
typedef long int delphi_integer;
typedef double   delphi_real;
#endif

#ifdef LD
typedef long int    delphi_integer;
typedef long double delphi_real;
#endif

/**
 * flag to debug the function(s) of reading SIZE file
 */
//#define DEBUG_IO_SIZE

#ifdef DEBUG_IO_SIZE
#define DEBUG_IO_FORCE
#endif

/**
 * flag to debug the function(s) of reading CHARGE file
 */
//#define DEBUG_IO_CHARGE

#ifdef DEBUG_IO_CHARGE
#define DEBUG_IO_FORCE
#endif

/**
 * flag to debug the function(s) of reading PDB file
 */
//#define DEBUG_IO_PDB

/**
 * flag to debug the function(s) of IDataMarshal class
 */
//#define DEBUG_DATAMARSHAL

/**
 * flag to debug the function(s) of CDelphiData class
 */
//#define DEBUG_DELPHI_MAP

/**
 * flag to debug the function(s) of CSpace class
 */
//#define DEBUG_DELPHI_SPACE

/**
 * flag to debug the function(s) of CDelphiFastSOR class
 */
//#define DEBUG_DELPHI_SOLVER

#ifdef DEBUG_DELPHI_SOLVER
#define DEBUG_DELPHI_SOLVER_MKDBSF1
#define DEBUG_DELPHI_SOLVER_MKDBSF
#define DEBUG_DELPHI_SOLVER_ITIT
#define DEBUG_DELPHI_SOLVER_RELFAC
#define DEBUG_DELPHI_SOLVER_SETBC
#define DEBUG_DELPHI_SOLVER_SETCRG
#endif

/**
 * flag to debug the function(s) of CDelphiEnergy class
 */
//#define DEBUG_DELPHI_ENERGY

/**
 * flag to debug the function(s) of CSite class
 */
//#define DEBUG_DELPHI_SITE

/**
 * flag to debug calling delphicpp in mcce
 */
//#define DEBUG_MCCE

/**
 * flag to debug constructing/destroying objects
 */
//#define DEBUG_OBJECT

}

using namespace delphi;

#endif //ENVIRONMENT_H_
