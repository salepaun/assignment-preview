#ifndef OPS_H
#define OPS_H

#include <Eigen/Dense>
#include <Eigen/Core>

#include "MathDefs.h"


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


inline static scalar calculateDistance(
    Vector2s const & _Xi, \
    Vector2s const & _Xj)
{
  scalar L = 0.0;
  Vector2s Xd = _Xi - _Xj;
  // Basic
  //L = sqrt((_Xi.x-_Xj.x)*(_Xi.x-_Xj.x) + (_Xi.y-_Xj.y)*(_Xi.y-_Xj.y));
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
  //L = (_Xi.x-_Xj.x)*(_Xi.x-_Xj.x) + (_Xi.y-_Xj.y)*(_Xi.y-_Xj.y);
  // More general
  //L = Xd.dot(Xd);
  // Using Eigen library
  L = Xd.squaredNorm();

  return L;
}

#endif // OPS_H
