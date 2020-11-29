#include "RigidBodyLCPCollisionResolver.h"

#include "../Ops.h"

using namespace std;


typedef set<RigidBodyCollision> T_Collisions;
typedef vector<RigidBody> T_RBCntr;




static bool processRB(int _Idx, RigidBody &_RB, MatrixXs &_A, VectorXs &_B)
{
  if (_RB.isFixed()) {
    return false;
  };

  return true;
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

  MatrixXs A(4,4);
  VectorXs b(4);

  T_Collisions::const_iterator I = rbcs.begin();
  for (int i=0; I != rbcs.end(); ++I, ++i) {
    b << I->nhat.dot(rbs[I->i0].computeWorldSpaceVelocityGivenPositionRelativeToCM(I->r0));
  }

#if MY_DEBUG > 0
  D1("B=" << b.transpose());
  D1("A=" << A);
#endif
}

