// SPDX-License-Identifier: LGPL-3.0-or-later
// Modified on branch nucleic-acids by noshou: SGrid operator+ defined inline so every file can inline it.

/**
 * @file misc_grid.h
 * @brief grid data types and associated operations
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef GRID_H_
#define GRID_H_

#include <iostream>

#include "../interface/environment.h"

using namespace std;
   

/**
 * New generic variable types for grid positions and value
 * SGrid <delphi_real>   : coord
 * SGrid <delphi_integer>: int_coord
 */
template <class T>  struct SGrid
{
   T nX,nY,nZ;
};

/**
 * SGridValue <delphi_real>    - grid_value
 * SGridValue <delphi_integer> - int_grid_value
 */
template <class T> struct SGridValue
{
   SGrid <T> nGrid;
   T nValue;
};   

/**
 * grid extreme values
 */
template <class T> struct SExtrema
{
   SGrid <T> nMin,nMax;
};

/**
 * double grid values
 */
struct SDoubleGridValue
{
   SGrid<delphi_real> fgCoord;
   delphi_real fVal1;
   delphi_real fVal2;
};

/*------------------- newly defined operators on grids ------------------*/

/*
 * overload operator <<
 */
template <class T> extern ostream& operator<<(ostream &, SGrid<T> &);
template <class T> extern ostream& operator<<(ostream &, SGridValue<T> &);
template <class T> extern ostream& operator<<(ostream &, SExtrema<T> &);

/*
 * overload operator +
 */
// Defined here, not in misc_grid_opts.cpp, so that every file can inline them (they run billions of times).
template <class T> inline SGrid <T> operator+(const SGrid <T> & nVector1, const SGrid <T> & nVector2)
{
   SGrid <T> nVector3;
   nVector3.nX = nVector1.nX + nVector2.nX;
   nVector3.nY = nVector1.nY + nVector2.nY;
   nVector3.nZ = nVector1.nZ + nVector2.nZ;
   return nVector3;
}
template <class T> inline SGrid <T> operator+(const SGrid <T> & nVector, const T & nValue)
{
   SGrid <T> nRtVector;
   nRtVector.nX = nVector.nX + nValue;
   nRtVector.nY = nVector.nY + nValue;
   nRtVector.nZ = nVector.nZ + nValue;
   return nRtVector;
}
template <class T> inline SGrid <T> operator+(const T & nValue, const SGrid <T> & nVector)
{
   SGrid <T> nRtVector;
   nRtVector.nX = nValue + nVector.nX;
   nRtVector.nY = nValue + nVector.nY;
   nRtVector.nZ = nValue + nVector.nZ;
   return nRtVector;
}

/*
 * overload operator -
 */
template <class T> extern SGrid <T> operator-(const SGrid <T> &, const SGrid <T> &);
template <class T> extern SGrid <T> operator-(const SGrid <T> &, const T &);
template <class T> extern SGrid <T> operator-(const T &, const SGrid <T> &);
template <class T> extern SGrid <T> operator-(const SGrid <T> &);
   
/*
 * overload operator *
 */
template <class T> extern SGrid <T> operator*(const SGrid <T> &, const SGrid <T> &);
template <class T> extern SGrid <T> operator*(const SGrid <T> &, const T &);
template <class T> extern SGrid <T> operator*(const T &, const SGrid <T> &);
   
/*
 * overload operator /
 */
template <class T, class N> extern SGrid <T> operator/(const SGrid <N> &, const T &);   
   
/*
 * new operator: dot product
 */
template <class T> extern T optDot(const SGrid <T> &, const SGrid <T> &);   
   
/*
 * new operator: cross product
 */
template <class T> extern SGrid <T> optCross(const SGrid <T> &, const SGrid <T> &); 
   
/*
 * new math func. optSQRT()
 */
extern SGrid <delphi_real> optSQRT(const SGrid <delphi_real> &);
   
/*
 * new math func. optABS()
 */
template <class T> extern SGrid <T> optABS(const SGrid <T> &);
         
/*
 * convert one type to another
 */
template <class T, class N> extern SGrid <T> optCast(const SGrid <N> &);
   
/*
 * round to the nearest delphi_integer
 */
extern SGrid <delphi_integer> optRound(const SGrid <delphi_real> &);
 
/*
 * new math func. optMin()
 */
template <class T> extern SGrid <T> optMin(const T &, const SGrid <T> &);
template <class T> extern SGrid <T> optMin(const SGrid <T> &, const T &);   
template <class T> extern SGrid <T> optMin(const SGrid <T> &, const SGrid <T> &);   
template <class T> extern T optMin(const SGrid <T> &);
   
/*
 * new math func. optMax()
 */
template <class T> extern SGrid <T> optMax(const T &, const SGrid <T> &);
template <class T> extern SGrid <T> optMax(const SGrid <T> &, const T &);
template <class T> extern SGrid <T> optMax(const SGrid <T> &, const SGrid <T> &);
template <class T> extern T optMax(const SGrid <T> &);

/*
 * new func. optMinSign()
 */
template <class T> extern SGrid <T> optMinSign(const T &, const SGrid <T> &);
template <class T> extern SGrid <T> optMinSign(const SGrid <T> &, const T &);

/*
 * new func. optMaxSign()
 */
template <class T> extern SGrid <T> optMaxSign(const T &, const SGrid <T> &);
template <class T> extern SGrid <T> optMaxSign(const SGrid <T> &, const T &);
  
/*
 * new func. optSubMin()
 */
template <class T> extern SGrid <T> optSubMin(const SGrid <T> &, const SGrid <T> &, const SGrid <T> &);
   
/*
 * new func. optSubMax()
 */
template <class T> extern SGrid <T> optSubMax(const SGrid <T> &, const SGrid <T> &, const SGrid <T> &);

/*
 * new func. optSum()
 */
template <class T> extern T optSum(const SGrid <T> &);

/*
 * new func. optComp()
 */
template <class T> extern T optComp(const SGrid <T> &, const int &);

/*
 * new func. optORLT()
 */
template <class T> extern bool optORLT(const SGrid <T> &, const SGrid <T> &);
template <class T> extern bool optORLT(const SGrid <T> &, const T &);
   
/*
 * new func. optANDLT()
 */
template <class T> extern bool optANDLT(const SGrid <T> &, const SGrid <T> &);  
template <class T> extern bool optANDLT(const SGrid <T> &, const T &);   
   
/*
 * new func. optORLE()
 */
template <class T> extern bool optORLE(const SGrid <T> &, const SGrid <T> &);
template <class T> extern bool optORLE(const SGrid <T> &, const T &);
  
/*
 * new func. optANDLE()
 */
template <class T> extern bool optANDLE(const SGrid <T> &, const SGrid <T> &);
template <class T> extern bool optANDLE(const SGrid <T> &, const T &);   
   
/*
 * new func. optORGT()
 */
template <class T> extern bool optORGT(const SGrid <T> &, const SGrid <T> &);   
template <class T> extern bool optORGT(const SGrid <T> &, const T &);   
   
/*
 * new func. optORGE()
 */
template <class T> extern bool optORGE(const SGrid <T> &, const SGrid <T> &);   
template <class T> extern bool optORGE(const SGrid <T> &, const T &);   
   
/*
 * new func. optANDGT()
 */
template <class T> extern bool optANDGT(const SGrid <T> &, const SGrid <T> &);
template <class T> extern bool optANDGT(const SGrid <T> &, const T &);   
   
/*
 * new func. optANDGE()
 */
template <class T> extern bool optANDGE(const SGrid <T> &, const SGrid <T> &); 
template <class T> extern bool optANDGE(const SGrid <T> &, const T &);   
   
/*
 * overlaod logical operator !=
 */
template <class T> extern bool operator!=(const SGrid <T> &, const T &);   
template <class T> extern bool operator!=(const SGrid <T> &, const SGrid <T> &);
    


#ifdef SP
/*
 * Single-precision build (delphi_real = float): the code mixes SGrid<float> with double literals and variables
 * (2.0*v, v+0.5, v/2., ...). These exact-match overloads convert the double to float so that such expressions compile
 * unchanged. They are not templates, so they win over the generic templates above and cannot be ambiguous.
 */
inline SGrid<float> operator+(const SGrid<float>& g, double v) { SGrid<float> r; r.nX = g.nX + (float)v; r.nY = g.nY + (float)v; r.nZ = g.nZ + (float)v; return r; }
inline SGrid<float> operator+(double v, const SGrid<float>& g) { return g + v; }
inline SGrid<float> operator-(const SGrid<float>& g, double v) { SGrid<float> r; r.nX = g.nX - (float)v; r.nY = g.nY - (float)v; r.nZ = g.nZ - (float)v; return r; }
inline SGrid<float> operator*(const SGrid<float>& g, double v) { SGrid<float> r; r.nX = g.nX * (float)v; r.nY = g.nY * (float)v; r.nZ = g.nZ * (float)v; return r; }
inline SGrid<float> operator*(double v, const SGrid<float>& g) { return g * v; }
inline SGrid<float> operator/(const SGrid<float>& g, double v) { SGrid<float> r; r.nX = g.nX / (float)v; r.nY = g.nY / (float)v; r.nZ = g.nZ / (float)v; return r; }
inline bool optORLT (const SGrid<float>& g, double v) { return optORLT (g, (float)v); }
inline bool optORGT (const SGrid<float>& g, double v) { return optORGT (g, (float)v); }
inline bool optORLE (const SGrid<float>& g, double v) { return optORLE (g, (float)v); }
inline bool optORGE (const SGrid<float>& g, double v) { return optORGE (g, (float)v); }
inline bool optANDLT(const SGrid<float>& g, double v) { return optANDLT(g, (float)v); }
inline bool optANDGT(const SGrid<float>& g, double v) { return optANDGT(g, (float)v); }
inline bool optANDLE(const SGrid<float>& g, double v) { return optANDLE(g, (float)v); }
inline bool optANDGE(const SGrid<float>& g, double v) { return optANDGE(g, (float)v); }
#endif
#endif // GRID_H_

