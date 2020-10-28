#ifndef OPS_H
#define OPS_H

#include <Eigen/Dense>
#include <Eigen/Core>

#include "MathDefs.h"

#include <stdio.h>


typedef std::pair<Vector2s, Vector2s> T_VecPair;


inline static Vector2s extractVectorIdx(const VectorXs& _X, int _Idx)
{
  return Vector2s(_X[_Idx], _X[_Idx+1]);
}

inline static Vector2s extractVector(const VectorXs& _X, int _N)
{
  int Idx = _N << 1;
  return extractVectorIdx(_X, Idx);
}

inline Matrix2s extractMassMxIdx(const VectorXs& _M, int _Idx)
{
  Matrix2s Mn; Mn << _M[_Idx], 0, 0, _M[_Idx+1];
  return Mn;
}

inline Matrix2s extractMassMx(const VectorXs& _M, int _N)
{
  int Idx = _N << 1;
  return extractMassMxIdx(_M, Idx);
}

inline static void addVectorDirectIdx(VectorXs& _X, int _Idx, Vector2s const & _Xn)
{
  _X[_Idx++] += _Xn.x();
  _X[_Idx] += _Xn.y();
}

inline static void addVectorDirect(VectorXs& _X, int _N, Vector2s const & _Xn)
{
  int Idx = _N << 1;
  addVectorDirectIdx(_X, Idx, _Xn);
}

inline static void addVector(VectorXs &_V, int _N, Vector2s const & _Vn)
{
  addVectorDirect(_V, _N, Vector2s(_Vn.dot(Vector2s(1,0)), _Vn.dot(Vector2s(0,1))));
}

inline static void updateVectorDirectIdx(VectorXs& _X, int _Idx, Vector2s const & _Xn)
{
  _X[_Idx++] = _Xn.x();
  _X[_Idx] = _Xn.y();
}

inline static void updateVectorDirect(VectorXs& _X, int _N, Vector2s const & _Xn)
{
  int Idx = _N << 1;
  updateVectorDirectIdx(_X, Idx, _Xn);
}

inline static void updateVector(VectorXs &_V, int _N, Vector2s const & _Vn)
{
  updateVectorDirect(_V, _N, Vector2s(_Vn.dot(Vector2s(1,0)), _Vn.dot(Vector2s(0,1))));
}


inline static void updateHessianWithJacobianBlock( \
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


inline static scalar calculateDistance(
    Vector2s const & _Xi, \
    Vector2s const & _Xj)
{
  scalar L = 0.0;
  Vector2s Xd = _Xi - _Xj;
  // Basic
  //L = sqrt((_Xi.x()-_Xj.x())*(_Xi.x()-_Xj.x()) + (_Xi.y()-_Xj.y())*(_Xi.y()-_Xj.y()));
  // More general
  //L = sqrt(Xd.dot(Xd));
  // Using Eigen library
  L = Xd.norm();

  return L;
}


inline static scalar calculateSquaredDistance(
    Vector2s const & _Xi, \
    Vector2s const & _Xj)
{
  scalar L = 0.0;
  Vector2s Xd = _Xi - _Xj;
  // Basic
  // L = (_Xi.x()-_Xj.x())*(_Xi.x()-_Xj.x()) + (_Xi.y()-_Xj.y())*(_Xi.y()-_Xj.y());
  // More general
  //L = Xd.dot(Xd);
  // Using Eigen library
  L = Xd.squaredNorm();

  return L;
}


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
      Vector2s X = extractVectorIdx(_X, i);
      Vector2s V = extractVectorIdx(_V, i);
      Vector2s M = extractVectorIdx(_M, i);
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



#endif // OPS_H
