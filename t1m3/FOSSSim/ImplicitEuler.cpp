#include "ImplicitEuler.h"
#include <Eigen/Dense>

#include "Ops.h"

#include <iostream>


#define MY_DEBUG
#undef MY_DEBUG

using namespace std;



/**
 * Local (static) helper functions.
 */

inline static void dVAll( \
    int _NDof, \
    int _NP, \
    VectorXs const &_Xn, \
    VectorXs const &_Vn, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_GradU, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dt)
{
  // Calculating Vn+1
  MatrixXs LS = _M - _dt*(_dt*_Hx + _Hv);
  VectorXs RS = -_dt * _GradU - _M * (_V - _Vn);
  VectorXs dV = LS.fullPivLu().solve(RS);
  _V += dV;
}


inline static void dXAll( \
    int _NDof, \
    int _NP, \
    VectorXs const &_Xn, \
    VectorXs const &_Vn, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M, \
    VectorXs const &_GradU, \
    MatrixXs const &_Hx, \
    MatrixXs const &_Hv, \
    scalar const &_dt)
{
  // Calculating Xn+1
  _X = _Xn + _dt * _V;
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
static scalar g_sTime = 0.0;

static scalar g_sNewtonError = 10.0e-9;




/**
 * Implicit Euler Newton step.
 */
static void calculateImplicitEulerStep(
    TwoDScene &_Scene, \
    scalar _dt, \
    int _NDof, \
    int _NP, \
    VectorXs const &_Xn, \
    VectorXs const &_Vn, \
    VectorXs &_Vi, \
    VectorXs &_X, \
    VectorXs &_V, \
    MatrixXs const &_M)
{
  VectorXs dX = _dt * _V; // !!!!! evaluation at (qn+h*vn+1_i, vn+1_i)
  VectorXs dV = VectorXs::Zero(_NDof);
  MatrixXs Hx = MatrixXs::Zero(_NDof, _NDof);
  MatrixXs Hv = MatrixXs::Zero(_NDof, _NDof);
  VectorXs GradU = VectorXs::Zero(_NDof);

  _Vi = _V;

  _Scene.accumulateGradU(GradU, dX, dV);
  _Scene.accumulateddUdxdx(Hx, dX, dV);
  _Scene.accumulateddUdxdv(Hv, dX, dV);

  handleFixedParticles(_Scene, _NDof, _NP, _X, _V, GradU, Hx, Hv);


#ifndef NDEBUG
#ifdef MY_DEBUG
  cout << "dX:" << dX << endl;
  cout << "dV:" << dV << endl;
  cout << "GradU:" << GradU << endl;
  cout << "Hx:" << Hx << endl;
  cout << "Hv:" << Hv << endl;
#endif
#endif


  // Calculating Vn+1
  dVAll(_NDof, _NP, _Xn, _Vn, _X, _V, _M, GradU, Hx, Hv, _dt);

  // Calculating Xn+1
  dXAll(_NDof, _NP, _Xn, _Vn, _X, _V, _M, GradU, Hx, Hv, _dt);
}




static bool keepNewton( \
    int _NDof, int _NP, \
    int &_Iter, \
    VectorXs const &_Vi, \
    VectorXs const &_Vi1)
{
  /*
  for(int i=0, n=0; i < _NDof; i+=2, n++) {
    scalar err = abs((_Vi.segment<2>(i) - _Vi1.segment<2>(i)).norm());
    if(err >= g_sNewtonError) {
#ifndef NDEBUG
      cout << __FUNCTION__ << "error[" << i << "]:" << err << endl;
      cout << "Vi+1 - Vi:  " << (_Vi1 - _Vi).transpose() << endl;
#endif
      return true;
    }
  }
  */

  scalar Err = (_Vi - _Vi1).squaredNorm();
#ifndef NDEBUG
  cout << (_Iter ? "  " : "**") \
    << (Err >= g_sNewtonError ? " keep " : " last ") \
    << "[" << _Iter << "] (Vi+1-Vi).norm()=" << Err << ", Vi+1 - Vi:  " << (_Vi1 - _Vi).transpose() << endl;
#endif
  ++_Iter;
  return Err >= g_sNewtonError;
}




/**
 * Main API
 */
bool ImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
  VectorXs& X = scene.getX();
  VectorXs& V = scene.getV();
  const VectorXs& M = scene.getM();
  assert(X.size() == V.size());
  assert(X.size() == M.size());

  // Implement implicit Euler here for extra credit!

  int NDof = X.size(), NP = scene.getNumParticles();
  assert(NDof%2 == 0);

  VectorXs Xn1 = X;
  VectorXs Vn1 = V;
  VectorXs Vi = V;

  // Mass matrix
  MatrixXs Mm = M.asDiagonal();

  FILE *Fd = energyDumpInit(g_sTime);

  int Iter = 0;
  do {
    calculateImplicitEulerStep( \
        scene, dt, \
        NDof, NP, \
        X, V, Vi, \
        Xn1, Vn1, Mm);
    
  } while (keepNewton(NDof, NP, Iter, Vi, Vn1));

  // Final assignments
  X = Xn1;
  V = Vn1;

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




static void sampleCode()
{
  /*
  /////////////
  // Some examples of working with vectors, matrices, etc.
  
  // How to get the force Jacobian from two d scene
  int ndof = x.size();
  assert( ndof%2 == 0 );
  // Note that the system's state is passed to two d scene as a change from the last timestep's solution
  VectorXs dx = dt*v;
  VectorXs dv = VectorXs::Zero(ndof);
  MatrixXs A = MatrixXs::Zero(ndof,ndof);
  scene.accumulateddUdxdx(A,dx,dv);
  scene.accumulateddUdxdv(A,dx,dv);

  // Some useful vector operations
  Vector2s egvec(1.0,0.0);
  // Outer product
  Matrix2s egmat = egvec*egvec.transpose();

  // Manipulating blocks of matrix
  MatrixXs egmat2(20,20);
  // 2x2 block that begins at position 3,4
  egmat2.block(3,4,2,2) += egmat;

  // Add a vector to a matrix's diagonal
  VectorXs egvec2(20);
  egvec2.setConstant(2.0);
  egmat2.diagonal() += egvec;

  // Invert a 10 x 10 matrix
  MatrixXs B(10,10);
  // Fill the matrix with random numbers
  B.setRandom();
  // Create a vector of length 10
  VectorXs b(10);
  // Fill the vector with random numbers
  b.setRandom();
  // Compute the solution to A*x = b
  VectorXs sln = B.fullPivLu().solve(b);
  // Verify that we computed the solution
  //std::cout << (B*sln-b).norm() << std::endl;  

  // End of example code.
  /////////////
  */
}
