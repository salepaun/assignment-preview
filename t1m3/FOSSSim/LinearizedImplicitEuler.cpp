#include "LinearizedImplicitEuler.h"
#include "Ops.h"

#include <iostream>


#define MY_DEBUG

using namespace std;



/**
 * Local (static) helper functions.
 */

inline static void dVAll( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_F, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  VectorXs V1 = _M.fullPivLu().solve(-_dT*_F);
  _V += V1;
}


inline static void dXAll( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_F, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  // Calculating Xn+1
  VectorXs X1 = _dT * _V;
  _X += X1;
}


inline static void dVStep( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_F, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  for(int i=0, n=0; n < _NP; n++) {
    i = n << 1;
    Vector2s V1n = _M.block<2,2>(i,i).fullPivLu().solve(-_dT*_F.segment<2>(i));
    _V.segment<2>(i) += V1n;
  };
}


inline static void dXStep( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_F, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  for(int i=0, n=0; n < _NP; n++) {
    i = n << 1;
    Vector2s X1n = _dT * _V.segment<2>(i);
    _X.segment<2>(i) += X1n;
  }
}




/**
 * Global variables.
 */
bool g_bUseStep = false;




/**
 * Main API
 */
bool LinearizedImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
  VectorXs& X = scene.getX();
  VectorXs& V = scene.getV();
  VectorXs const & M = scene.getM();
  assert(X.size() == V.size());
  assert(X.size() == M.size());

  // Your code goes here!
  int NDof = X.size(), NP = scene.getNumParticles();
  assert(NDof%2 == 0);
  VectorXs dX = dt * V;
  VectorXs dV = VectorXs::Zero(NDof);
  MatrixXs Hx = MatrixXs::Zero(NDof, NDof);
  MatrixXs Hv = MatrixXs::Zero(NDof, NDof);
  VectorXs F = VectorXs::Zero(NDof);

#ifndef NDEBUG
#ifdef MY_DEBUG
  cout << "Calling Udxdx" << endl;
  scene.accumulateddUdxdx(Hx, dX, dV);
  cout << "Calling Udxdv" << endl;
  scene.accumulateddUdxdv(Hv, dX, dV);
  cout << "Calling GradU" << endl;
  scene.accumulateGradU(F, dX, dV);

  cout << "dX:" << dX << endl;
  cout << "dV:" << dV << endl;
  cout << "F:" << F << endl;
  cout << "Hx:" << Hx << endl;
  cout << "Hv:" << Hv << endl;
#endif
#endif

  // 
  // Mass matrix
  MatrixXs M1 = MatrixXs::Zero(NDof, NDof);
  for(int i=0, j=0; i < NDof; ++i, ++j)
    M1(i,j) = M(i);

  M1 -= dt*(dt*Hx + Hv);
  // tests
  // MatrixXs MI = M1.inverse();
  // VectorXs V1 = -dt * MI * F;

  if(g_bUseStep && NP > 100) {
    // Works only when particles don't have common Hessian (each has it's own forces),
    // which was the case for all previous assignments
 
    // Calculating Vn+1
    dVStep(NDof, NP, X, V, M1, F, Hx, Hv, dt);
    // Calculating Xn+1
    dXStep(NDof, NP, X, V, M1, F, Hx, Hv, dt);
  } else {
    // Might require more CPU for more particles (??)
    // but can calculate relations between particles
 
    // Calculating Vn+1
    dVAll(NDof, NP, X, V, M1, F, Hx, Hv, dt);
    // Calculating Xn+1
    dXAll(NDof, NP, X, V, M1, F, Hx, Hv, dt);
  };

#ifndef NDEBUG
#ifdef MY_DEBUG
  dumpParticles(X, V, M, 0, 0, true, true);
#endif
#endif
  
  return true;
}
