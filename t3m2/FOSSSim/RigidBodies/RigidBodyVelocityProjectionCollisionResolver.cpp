#include "RigidBodyVelocityProjectionCollisionResolver.h"

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

#include "../Ops.h"





typedef std::set<RigidBodyCollision> T_Collisions;
typedef std::vector<RigidBody> T_RBCntr;
typedef std::map<int, T_RBCntr::iterator> T_FreeRBs;




static int g_DebugCntrLimit = 10;



static std::ostream & operator << (std::ostream &_s, T_FreeRBs::value_type const &_o) {
  return _s << _o.first;
};






/**
 * Bulds Gamma matrix per contact point.
 * @param _K number of not fixed rigid bodies
 * @param _I index of this contact point first RB
 * @param _J index of this contact point second RB
 */
static void buildGamma(int _K,
    int _Ii, int _Ij,
    T_FreeRBs const &_FreeRBs,
    Vector2s const &_Ri,
    Vector2s const &_Rj,
    MatrixXs &_Gamma)
{
  T_FreeRBs::const_iterator RBiI, RBjI;

  _Gamma.resize(2, 3*_K);
  _Gamma.setZero();

  RBiI = _FreeRBs.find(_Ii);
  RBjI = _FreeRBs.find(_Ij);
  scalar ThetaI = RBiI != _FreeRBs.end() ? RBiI->second->getTheta() : 1.0;
  scalar ThetaJ = RBjI != _FreeRBs.end() ? RBjI->second->getTheta() : 0.0;

  T_FreeRBs::const_iterator I = _FreeRBs.begin();
  for (int i=0; I != _FreeRBs.end(); ++I, i+=3) {
    if (I->first == _Ii) {
      _Gamma.block<2,3>(0, i) << 
        1, 0, -sin(ThetaI)*_Ri.x() - cos(ThetaI)*_Ri.y(),
        0, 1,  cos(ThetaI)*_Ri.x() - sin(ThetaI)*_Ri.y();
    } else if (I->first == _Ij) {
      _Gamma.block<2,3>( 0, i) << 
        -1, 0, sin(ThetaJ)*_Rj.x() + cos(ThetaJ)*_Rj.y(),
        0, -1,-cos(ThetaJ)*_Rj.x() + sin(ThetaJ)*_Rj.y();
    };
  };

#if MY_DEBUG > 1
  D1(" ** ThetaI=" << ThetaI << ", ThetaJ=" << ThetaJ);
  D1(" ** Gamma=\n" << _Gamma);
#endif
}




/**
 * Builds matrix N for all collision points.
 */
static void buildN(
    int _K,
    T_FreeRBs const &_FreeRBs,
    T_Collisions const &_Colls,
    MatrixXs &_N)
{
  _N.resize(3*_K, _Colls.size());
  _N.setZero();

  T_Collisions::const_iterator I = _Colls.begin();
  for (int i=0; I != _Colls.end(); ++I, ++i) {
      int Ii = I->i0;
      int Ij = I->i1;
      MatrixXs Gamma;
      buildGamma(_K, Ii, Ij, _FreeRBs, I->r0, I->r1, Gamma);
      _N.col(i) = Gamma.transpose() * I->nhat;
  }

#if MY_DEBUG > 1
  D1(" ** N=\n" << _N);
#endif
}





/**
 * Builds M matrix.
 */
static void buildM(
    int _K,
    T_FreeRBs const &_FreeRBs,
    MatrixXs &_M)
{
  _M.resize(3*_K, 3*_K);
  _M.setZero();

  T_FreeRBs::const_iterator I = _FreeRBs.begin();
  for (int i=0; I != _FreeRBs.end(); ++I, i+=3) {
    _M.block<3,3>(i,i).diagonal() << I->second->getM(), I->second->getM(), I->second->getI();
  };

#if MY_DEBUG > 1
  D1(" ** M=\n" << _M);
#endif
}




/**
 * Builds Qdot vector.
 */
static void buildQdot(
    int _K,
    T_FreeRBs const &_FreeRBs,
    VectorXs &_Qdot)
{
  _Qdot.resize(3*_K);
  _Qdot.setZero();

  T_FreeRBs::const_iterator I = _FreeRBs.begin();
  for (int j=0; I != _FreeRBs.end(); ++I, j+=3) {
    _Qdot.segment<3>(j) << I->second->getV()(0), I->second->getV()(1), I->second->getOmega();
  }

#if MY_DEBUG > 1
  D1(" ** Qdot=\n" << _Qdot.transpose());
#endif
}




/**
 * Builds container with free RigidBodies.
 */
static void buildFreeRBs(
    T_RBCntr &_RBs,
    T_FreeRBs &_FreeRBs,
    int &_K)
{
  _FreeRBs.clear();

  T_RBCntr::iterator I = _RBs.begin();
  for (int i=0; I != _RBs.end(); ++I, i++) {
    if (!I->isFixed()) {
      _FreeRBs.insert(T_FreeRBs::value_type(i, I));
    };
  };

  _K = _FreeRBs.size();

#if MY_DEBUG > 1
  dumpContainer(g_DebugCntrLimit, std::cout,
      __FUNCTION__, ": FreeRBs:",
      _FreeRBs.size(), _FreeRBs.begin(), _FreeRBs.end()) << std::endl;
#endif
}




void RigidBodyVelocityProjectionCollisionResolver::resolveCollisions( std::vector<RigidBody>& rbs, const std::set<RigidBodyCollision>& rbcs )
{
  // Your code goes here!

  /*
  // Example of using QuadProg++
  // For detailed documentation, please see FOSSSim/quadprog/QuadProg++.hh
 
  // Matrix in quadratic form of objective
  QuadProgPP::Matrix<scalar> G;
  // Equality constraints
  QuadProgPP::Matrix<scalar> CE;
  // Inequality constraints
  QuadProgPP::Matrix<scalar> CI;
  
  QuadProgPP::Vector<scalar> g0;
  QuadProgPP::Vector<scalar> ce0;
  QuadProgPP::Vector<scalar> ci0;
  QuadProgPP::Vector<scalar> x;

  // M = 16  0 0
  //      0 16 0
  //      0  0 289.75
  G.resize(3,3);
  for(int i=0; i< 3; i++) for(int j=0; j< 3; j++) G[i][j]=0;
  G[0][0] = 16; 
  G[1][1] = 16;
  G[2][2] = 289.75;

  // -M \dot q = -0
  //             139.52
  //            -0
  g0.resize(3);
  g0[0] = 0;
  g0[1] = 139.52;
  g0[2] = 0;

  // No equality constraints, currently
  CE.resize(3,0);
  ce0.resize(0);

  // Compute the number of inequality constraints
  CI.resize(3,24);

  MatrixXs tempN(3,24);
  tempN << 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           1,  1,  1,  1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          -5, -4, -3, -2, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 24; ++j ) CI[i][j] = tempN(i,j);

  // Constant term added to inequality constraints
  ci0.resize(24);
  for(int i=0; i< 24; i++) ci0[i] = 0;

  // Solution
  x.resize(3);

  solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
  */

  if (rbcs.size())
  {
    MatrixXs M;
    MatrixXs N;
    VectorXs Qdot;
    T_FreeRBs FreeRBs;
    int K = 0;

    buildFreeRBs(rbs, FreeRBs, K);
    buildQdot(K, FreeRBs, Qdot); 
    buildM(K, FreeRBs, M);
    buildN(K, FreeRBs, rbcs, N);

    VectorXs b(M.rows());

    b = (-M * Qdot).transpose();

    QuadProgPP::Matrix<scalar> G;
    // Equality constraints
    QuadProgPP::Matrix<scalar> CE;
    // Inequality constraints
    QuadProgPP::Matrix<scalar> CI;

    QuadProgPP::Vector<scalar> g0;
    QuadProgPP::Vector<scalar> ce0;
    QuadProgPP::Vector<scalar> ci0;
    QuadProgPP::Vector<scalar> x;

    int D = 3*K;

    G.resize(D,D);
    g0.resize(D);

    CE.resize(D,0);
    ce0.resize(0);

    CI.resize(D, N.cols());
    ci0.resize(N.cols());

    x.resize(D);

    for(int i=0; i < M.rows(); ++i) G[i][i] = M(i,i);
    for(int i=0; i < b.size(); ++i) g0[i] = b[i];

    for(int i=0; i < N.cols(); ++i) ci0[i] = 0;
    for(int i=0; i < D; ++i) for(int j=0; j < N.cols(); ++j) CI[i][j] = N(i,j);


#if MY_DEBUG > 0
    if (rbcs.size()) {
      D1("** Collisions=" << rbcs.size()
          << "\n** g0=" << b.transpose()
          << "\n** G=\n" << M
          << "\n** ci0=" << ci0
          << "\n** CI=\n" << CI);
    }
#endif

    double res = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);


#if MY_DEBUG > 0
    if (rbcs.size()) {
      D1("** Collisions=" << rbcs.size()
          << ", solve_quadprog res=" << res
          << "\n** v=\n" << x);
    }
#endif

    if (!isinf(res)) {
      T_FreeRBs::iterator I = FreeRBs.begin();
      for (int i=0; I != FreeRBs.end(); ++I, i+=3) {
        I->second->getV()[0] = x[i];
        I->second->getV()[1] = x[i+1];
        I->second->getOmega() = x[i+2];
      };
    };
  };
}
