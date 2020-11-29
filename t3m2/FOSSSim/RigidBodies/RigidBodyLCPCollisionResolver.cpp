#include "RigidBodyLCPCollisionResolver.h"

#include "../Ops.h"



using namespace std;


typedef set<RigidBodyCollision> T_Collisions;
typedef vector<RigidBody> T_RBCntr;






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

  int D = rbcs.size();
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
    J = I;
    for (int j=0; J != rbcs.end(); ++J, ++j) {
      A(i,j) = I->nhat.dot(J->nhat) / rbs[I->i0].getM() + I->nhat.dot(J->r0) * cross2s(J->nhat, (I->r0 / rbs[I->i0].getI()));
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
  dOmegas.resize(rbs.size());
  I = rbcs.begin();
  for (int i=0; I != rbcs.end(); ++I, ++i) {
    if (lambda(i) > 0) {
      int i0 = I->i0;
      int i1 = I->i1;
      RigidBody &A = rbs[i0];
      RigidBody &B = rbs[i1];
      Vector2s par1 = -2 *lambda(i) * I->nhat;
      dVs[i0] += par1;
      dVs[i1] -= par1;
    };
  };

  vector<Vector2s>::iterator dVsIter = dVs.begin();
  for (int i=0; dVsIter != dVs.end(); ++dVsIter, ++i) {
    rbs[i].getV() += *dVsIter / rbs[i].getM();
  };
}

