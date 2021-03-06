#include "LinearizedImplicitEuler.h"
#include "Ops.h"

#include <iostream>


#define MY_DEBUG
//#undef MY_DEBUG

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
    VectorXs const &_GradU, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  VectorXs V1 = _M.fullPivLu().solve(-_dT*_GradU);
  _V += V1;
}


inline static void dXAll( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_GradU, \
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
    VectorXs const &_GradU, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dT)
{
  for(int i=0, n=0; n < _NP; n++) {
    i = n << 1;
    Vector2s V1n = _M.block<2,2>(i,i).fullPivLu().solve(-_dT*_GradU.segment<2>(i));
    _V.segment<2>(i) += V1n;
  };
}


inline static void dXStep( \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_GradU, \
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


inline static void handleFixedParticles( \
    TwoDScene const &_Scene, \
    int _NDof, \
    int _NP, \
    VectorXs &_X, \
    VectorXs &_V, \
    VectorXs &_GradU, \
    MatrixXs &_Hx, \
    MatrixXs &_Hv)
{
  for(int n=0, i=0, j=0; n < _NP; ++n) {
    if(_Scene.isFixed(n)) {
      i = n << 1;
      j = i + 1;
      _V.segment<2>(i).setZero();
      _GradU.segment<2>(i).setZero();
      _Hx.col(i).setZero(); _Hx.row(i).setZero(); _Hx(i,i) = 1;
      _Hx.col(j).setZero(); _Hx.row(j).setZero(); _Hx(j,j) = 1;
      _Hv.col(i).setZero(); _Hv.row(i).setZero(); _Hv(i,i) = 1;
      _Hv.col(j).setZero(); _Hv.row(j).setZero(); _Hv(j,j) = 1;
    };
  };
}




/**
 * Global variables.
 */
static bool g_bUseStep = false;
static scalar g_sTime = 0.0;




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
  VectorXs GradU = VectorXs::Zero(NDof);

  scene.accumulateGradU(GradU, dX, dV);
  scene.accumulateddUdxdx(Hx, dX, dV);
  scene.accumulateddUdxdv(Hv, dX, dV);

  handleFixedParticles(scene, NDof, NP, X, V, GradU, Hx, Hv);


#ifndef NDEBUG
#ifdef MY_DEBUG
  cout << "dX:" << dX << endl;
  cout << "dV:" << dV << endl;
  cout << "GradU:" << GradU << endl;
  cout << "Hx:" << Hx << endl;
  cout << "Hv:" << Hv << endl;
#endif
#endif

  FILE *Fd = energyDumpInit(g_sTime);

  // 
  // Mass matrix
  MatrixXs M1 = MatrixXs::Zero(NDof, NDof);
  for(int i=0; i < NDof; ++i)
    M1(i,i) = M(i);

  M1 -= dt*(dt*Hx + Hv);
  // tests
  // MatrixXs MI = M1.inverse();
  // VectorXs V1 = dt * MI * GradU;

  if(g_bUseStep && NP > 100) {
    // Works only when particles don't have common Hessian (each has it's own forces),
    // which was the case for all previous assignments
 
    // Calculating Vn+1
    dVStep(NDof, NP, X, V, M1, GradU, Hx, Hv, dt);
    // Calculating Xn+1
    dXStep(NDof, NP, X, V, M1, GradU, Hx, Hv, dt);
  } else {
    // Might require more CPU for more particles (??)
    // but can calculate relations between particles
 
    // Calculating Vn+1
    dVAll(NDof, NP, X, V, M1, GradU, Hx, Hv, dt);
    // Calculating Xn+1
    dXAll(NDof, NP, X, V, M1, GradU, Hx, Hv, dt);
  };


  //
  // Scene kinetic energy
  scalar E = scene.computeKineticEnergy();
  Vector2s VTotal = Vector2s::Zero();
  for(int i=0; i < NDof; i+=2) {
    VTotal += V.segment<2>(i);
  };
  energyDumpEnd(g_sTime, dt, Fd, E, VTotal);


#ifndef NDEBUG
  dumpParticles(X, V, M, 0, 0, true, true);
#endif
  
  return true;
}
