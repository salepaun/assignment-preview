#include "RigidBodyLCPCollisionResolver.h"

#include <vector>
#include <map>
#include <algorithm>

#include "../Ops.h"


using namespace std;


typedef set<RigidBodyCollision> T_Collisions;
typedef vector<RigidBody> T_RBCntr;
typedef map<int, T_RBCntr::iterator> T_FreeRBs;




int g_DebugCntrLimit = 10;



ostream & operator << (ostream &_s, T_FreeRBs::value_type const &_o) {
  return _s << _o.first;
};






/**
 * Bulds Gamma matrix per contact point.
 * @param _K number of not fixed rigid bodies
 * @param _I index of this contact point first RB
 * @param _J index of this contact point second RB
 */
void buildGamma(int _K,
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
  scalar ThetaI = RBiI != _FreeRBs.end() ? RBiI->second->getTheta() : 0.0;
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
        0, -1,-cos(ThetaJ)*_Rj.x() - sin(ThetaJ)*_Rj.y();
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
void buildN(
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
void buildM(
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
void buildQdot(
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
void buildFreeRBs(
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
  dumpContainer(g_DebugCntrLimit, cout,
      __FUNCTION__, ": FreeRBs:",
      _FreeRBs.size(), _FreeRBs.begin(), _FreeRBs.end()) << endl;
#endif
}






void RigidBodyLCPCollisionResolver::resolveCollisions( std::vector<RigidBody>& rbs, const std::set<RigidBodyCollision>& rbcs )
{
  // Your code goes here!

  // Example of using the lcp solver
  /*
  MatrixXs A(4,4);
  A << 0.828869, 0.337798, -0.28125, -0.21875,
       0.337798, 0.828869, -0.21875, -0.28125,
       -0.28125, -0.21875, 0.828869, 0.337798,
       -0.21875, -0.28125, 0.337798, 0.828869;

  VectorXs b(4);
  b << -1, -1, -1, -1;

  VectorXs lambda(4);
  lcputils::solveLCPwithODE( A, b, lambda );
  */

  /*
  int D = rbcs.size() << 1;
  MatrixXs A(D, D);
  VectorXs b(D);
  VectorXs lambda(D);
  vector<Vector2s> dVs;
  vector<scalar> dOmegas;

  T_Collisions::const_iterator I,J;


  I = rbcs.begin();
  for (int i=0; I != rbcs.end(); ++I, ++i) {
    b << I->nhat.dot(rbs[I->i0].computeWorldSpaceVelocityGivenPositionRelativeToCM(I->r0));
  }

  I = rbcs.begin();
  for (int i=0; I != rbcs.end(); ++I, ++i) {
    J = rbcs.begin();
    for (int j=0; J != rbcs.end(); ++J, ++j) {
      int i0 = I->i0;
      int i1 = I->i1;
      RigidBody &Arb = rbs[i0];
      RigidBody &Brb = rbs[i1];
      A(i,j) = I->nhat.dot(J->nhat / Arb.getM()) + I->nhat.dot(J->r0) * cross2s(J->nhat, (I->r0 / Arb.getI()));
    };
  };

  lcputils::solveLCPwithODE(A, b, lambda);

#if MY_DEBUG > 0
  if (rbcs.size()) {
    D1("** Collisions=" << rbcs.size()
        << "\n** B=" << b.transpose()
        << "\n** L=" << lambda.transpose()
        << "\n** A=\n" << A);
  }
#endif

  dVs.resize(rbs.size());
  fill(dVs.begin(), dVs.end(), Vector2s::Zero());
  dOmegas.resize(rbs.size());
  fill(dOmegas.begin(), dOmegas.end(), 0.0);

  I = rbcs.begin();
  for (int i=0; I != rbcs.end(); ++I, ++i) {
    if (lambda(i) > 0) {
      int i0 = I->i0;
      int i1 = I->i1;
      RigidBody &A = rbs[i0];
      RigidBody &B = rbs[i1];
      Vector2s par1 = lambda(i) * I->nhat;
      dVs[i0] += par1;
      dVs[i1] -= par1;
    };
  };

  vector<Vector2s>::iterator dVsIter = dVs.begin();
  for (int i=0; dVsIter != dVs.end(); ++dVsIter, ++i) {
    rbs[i].getV() += *dVsIter / rbs[i].getM();
  };
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

    MatrixXs A(N.cols(), N.cols());
    VectorXs b(N.cols());
    VectorXs lambda(N.cols());
    vector<Vector2s> dVs;
    vector<scalar> dOmegas;

    MatrixXs Minv = M.inverse();

    A = N.transpose() * Minv * N;
    b = N.transpose() * Qdot;

    lcputils::solveLCPwithODE( A, b, lambda );

    VectorXs dV = Minv * N * lambda;

#if MY_DEBUG > 0
    if (rbcs.size()) {
      D1("** Collisions=" << rbcs.size()
          << "\n** B=" << b.transpose()
          << "\n** L=" << lambda.transpose()
          << "\n** dV=" << dV.transpose()
          << "\n** A=\n" << A);
    }
#endif

    T_FreeRBs::iterator I = FreeRBs.begin();
    for (int i=0; I != FreeRBs.end(); ++I, i+=3) {
      I->second->getV() += dV.segment<2>(i);
      I->second->getOmega() += dV(i+2);
    };
  };
}

