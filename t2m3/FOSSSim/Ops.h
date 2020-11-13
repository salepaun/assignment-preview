#ifndef OPS_H
#define OPS_H

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <set>
#include <utility>
#include <iostream>


#include <Eigen/Dense>

#include "MathDefs.h"

using namespace std;


/**
 * Type definitions.
 */

typedef pair<Vector2s, Vector2s> T_VecPair;



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



/**********************************************************************************
 * Local stream helper functions
 */


template<typename T_Iter>
ostream & dumpContainer(int _Limit, ostream &_s,
    char const * const _Function,
    char const * const _Name,
    size_t _Size,
    T_Iter const &_Begin,
    T_Iter const &_End)
{
  if (_Function) _s.write(_Function, strlen(_Function));
  if (_Name) _s.write(_Name, strlen(_Name));
  _s << "Size=" << _Size;
  T_Iter I = _Begin;
  for(int i=0; i<_Limit && I!=_End; ++I, ++i) {
    _s << ", " << *I;
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
static void combination2(int _nA, int _nB, set<pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < _nA; ++i)
    for(int j=0; j < _nB; ++j)
      _Set.insert(pair<int,int>(i,j));
}




/**
 * Creates combination '2' of n elements.
 * Clears the set.
 */
static void combination(int n, set<pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < n; ++i)
    for(int j=i+1; j < n; ++j)
      _Set.insert(pair<int,int>(i,j));
}




/**
 * Creates combination '2' of n elements.
 * Clears the set.
 */
static void combinationAll(int n, set<pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < n; ++i)
    for(int j=i; j < n; ++j)
      _Set.insert(pair<int,int>(i,j));
}



#endif // OPS_H
