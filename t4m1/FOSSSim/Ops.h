#ifndef OPS_H
#define OPS_H

#include <stdio.h>
#include <string.h>
#include <set>
#include <utility>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>


#include <Eigen/Dense>

#include "MathDefs.h"


#ifndef NDEBUG
#include <assert.h>

// For now: 0,1,2,3,4,5
#define MY_DEBUG 3
// 0,1,2
#define MY_TIMING 1


#else
#undef MY_DEBUG
#undef MY_TIMING
#endif


#ifndef NDEBUG

#define D1(a) std::cout << __FUNCTION__ << ":" << a << std::endl;

#define D_TIME_SEC(a,b) std::cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< std::endl; a = clock();
#define D_TIME_SEC2(a,b,c) std::cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< std::endl; a = clock();
#define D_TIME_SEC3(a,b,c,d) std::cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< ": " << d \
<< std::endl; a = clock();
#define D_TIME_SEC4(a,b,c,d,e) std::cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< ": " << d \
<< ": " << e \
<< std::endl; a = clock();

#endif


#define MIN(a,b) a<b ? a : b
#define MAX(a,b) a>b ? a : b




/**
 * Type definitions.
 */

typedef std::pair<Vector2s, Vector2s> T_VecPair;



/**********************************************************************************
 * Local CGI related helper functions
 */


static void dumpParticles( \
    VectorXs const & _X, \
    VectorXs const & _V, \
    VectorXs const & _M, \
    int _Max=0, \
    int _Min=0, \
    bool _WithV=false, \
    bool _Append=false, \
    char const * const _FileName="ParticlesDump.xml")
{
  FILE *fd = NULL;
  fd = _Append ? fopen(_FileName, "a") : fopen(_FileName, "w"); 
  if(!fd) {
    perror("Failed to open scene dump file");
  } else {
    int Size = _Max==0 ? _X.size() : (_Max << 1) + 1;
    for(int i=0, n=0; i < Size; i+=2, ++n) {
      Vector2s X = _X.segment<2>(i);
      Vector2s V = _V.segment<2>(i);
      Vector2s M = _M.segment<2>(i);
      if(!_WithV)
        V.setZero();
      fprintf(fd, \
          "<particle m=\"%.4lf\" px=\"%.4lf\" py=\"%.4lf\" vx=\"%.4lf\" vy=\"%.4lf\" fixed=\"0\"/>\n", \
          M.x(), X.x(), X.y(), V.x(), V.y());
    };
    fclose(fd);
  };
}




static FILE * energyDumpInit(scalar const &_sTime)
{
  FILE *fd = NULL;

  if(_sTime) {
    fd = fopen("./kinetic_energy_output.txt", "a");
  } else {
    fd = fopen("./kinetic_energy_output.txt", "w");
    if(fd)
      fprintf(fd, "# Time\t Kinetic Energy\t Vx\t Vy\n");
    else
      perror("Failed to open kinetic energy dump");
  };

  return fd;
}



static void energyDumpEnd(scalar &_sTime, scalar const &_Dt, FILE *_Fd, scalar const &_E, Vector2s const &_V)
{
  _sTime += _Dt;

  if(_Fd) {
    fprintf(_Fd, "%'.2lf\t%'.5lf\t%'.5lf\t%'.5lf\n", _sTime, _E, _V[0], _V[1]);
    fclose(_Fd);
  };
}



static void findAlpha( \
    Vector2s const &_X1, \
    Vector2s const &_X2, \
    Vector2s const &_X3, 
    Vector2s const &_V1, \
    Vector2s const &_V2, \
    Vector2s const &_V3, 
    Vector2s &_Ne, \
    Vector2s &_Ve, \
    scalar &_sAlpha, \
    bool _bDebug=true)
{
    Vector2s A = _X1 - _X2;
    Vector2s B = _X3 - _X2;
    Vector2s EAxN = B;
    Vector2s EPAx; // Edge Particle axis
    EAxN.normalize();

    // alpha parameter (edge drop point absolute distance)
    scalar B2Ln = B.squaredNorm();
    _sAlpha = A.dot(B) / B2Ln;

    // 
    // Alpha constraints
    //
    if (_sAlpha <= 0) {
      // away from or on X2 edge vertex
      _sAlpha = 0;
      EPAx = (_X2 - _X1);
      _Ve = _V2;
    } else if (_sAlpha >= 1) {
      // away from or on X3 edge vertex
      _sAlpha = 1;
      EPAx = (_X3 - _X1);
      _Ve = _V3;
    } else {
      // calculating Normal and distance with 'drop' point
      EPAx = _sAlpha * B - A;
      _Ve = _V2 + _sAlpha * (_V3 - _V2);
    };

    _Ne = EPAx;

#if MY_DEBUG > 4
    if (_bDebug) {
      D1("(): Alpha:" << _sAlpha \
        << ", A:(" << A.transpose() << "), B:(" << B.transpose() << "), N:(" << _Ne.transpose() << ")");
    };
#endif
}




inline static void updateHessianWithBlock( \
    int _Ii, \
    int _Ij, \
    MatrixXs &_H, \
    Matrix2s const &_K)
{
  _H.block<2,2>(_Ii,_Ii) += _K;
  _H.block<2,2>(_Ii,_Ij) += -_K;
  _H.block<2,2>(_Ij,_Ii) += -_K;
  _H.block<2,2>(_Ij,_Ij) += _K;
}







/**********************************************************************************
 * Local container helper functions
 */

/**
 * Finds 2 vectors insersection.
 * Vectors have to be sorted
 * @param _A a smaller sized vector
 * @param _B a bigger sized vector
 * @param _C destination
 */
template<typename T>
inline bool findSortedCntrIntersect(T const &_A, T const &_B, T &_C) {
  typename T::const_iterator IA;
  typename T::const_iterator IB1, IB2;
  register int OrigSizeC = _C.size();
  for (IA=_A.cbegin(), IB1=_B.cbegin(); IA!=_A.cend(); ++IA) {
    IB2 = find(IB1, _B.cend(), (*IA));
    if (IB2 != _B.cend()) {
      _C.push_back((*IA));
      IB1 = IB2;
    };
  };

  return _C.size() - OrigSizeC;
}


/**
 * Finds 2 vectors insersection.
 * Vectors will be sorted
 * @param _A a vector
 * @param _B a vector
 * @param _C destination
 */
template<typename T>
inline bool findCntrIntersect(T &_A, T &_B, T &_C) {
  register int SizeA, SizeB;
  SizeA = _A.size(); SizeB = _B.size();
  if(SizeA && SizeB) {
    sort(_A.begin(), _A.end());
    sort(_B.begin(), _B.end());

    if (SizeA < SizeB) {
      _C.reserve(SizeA);
      return findSortedCntrIntersect(_A, _B, _C);
    } else {
      _C.reserve(SizeB);
      return findSortedCntrIntersect(_B, _A, _C);
    };
  };

  return false;
}


/**********************************************************************************
 * Local stream helper functions
 */


  template<typename T_Iter>
std::ostream & dumpContainer(int _Limit, std::ostream &_s,
    char const * const _Function,
    char const * const _Name,
    size_t _Size,
    T_Iter const &_Begin,
    T_Iter const &_End,
    char const * const _Sep=", ")
{
  if (_Function) _s.write(_Function, strlen(_Function));
  if (_Name) _s.write(_Name, strlen(_Name));
  if (!_Limit) _Limit = _Size;
  _s << "Size=" << _Size;
  T_Iter I = _Begin;
  for(int i=0; i<_Limit && I!=_End; ++I, ++i) {
    _s << _Sep << *I;
  };

  return _s;
}




/**********************************************************************************
 * Local mathematical functions.
 */


/**
 * Returns factorial of n.
 */
unsigned int inline factorial(int n) {
  int F = 1;
  for(int i=2; i<=n; ++i) {
    F *= i;
  };

  return F;
}




/**
 * Calculates combinations k of n.
 */
unsigned int inline nCk(int n, int k)
{
  return factorial(n)/factorial(k)/factorial(n-k);
}




/**
 * Creates combination of '2' elements out of 2 sets.
 * Clears the set.
 *
 * @param _nA - the size of set A
 * @param _nB - the size of set B
 */
static void combination2(int _nA, int _nB, std::set<std::pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < _nA; ++i)
    for(int j=0; j < _nB; ++j)
      _Set.insert(std::pair<int,int>(i,j));
}




/**
 * Creates combination '2' of n elements.
 * Clears the set.
 */
static void combination(int n, std::set<std::pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < n; ++i)
    for(int j=i+1; j < n; ++j)
      _Set.insert(std::pair<int,int>(i,j));
}




/**
 * Creates combination '2' of n elements.
 * Clears the set.
 */
static void combinationAll(int n, std::set<std::pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < n; ++i)
    for(int j=i; j < n; ++j)
      _Set.insert(std::pair<int,int>(i,j));
}




/**********************************************************************************
 * Local eigen helper functions
 */


static inline scalar cross2s(Vector2s const &_A, Vector2s const &_B)
{
  return _A.x()*_B.y() - _B.x()*_A.y();
}



static Vector2s rotate(Vector2s const &_V, scalar _Theta)
{
  Matrix2s R; R << cos(_Theta), -sin(_Theta), sin(_Theta), cos(_Theta);
  return R * _V;
}



#endif // OPS_H
