#include "RigidBodyGravityForce.h"

#include "../Ops.h"


static Vector2s calculateSimpleGravity(T_RBodyCntr::const_iterator const &_I, Vector2s const &_g)
{
  return _I->getM() * _g;
}



scalar RigidBodyGravityForce::computePotentialEnergy( const std::vector<RigidBody>& rbs )
{
  // Your code goes here!
  scalar PE = 0.0;

  T_RBodyCntr::const_iterator I = rbs.begin();
  for (; I != rbs.end(); ++I) {
    PE += calculateSimpleGravity(I, m_g).dot(I->getX());
  };

  return -PE;
}



void RigidBodyGravityForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  // Your code goes here!
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...

  T_RBodyCntr::iterator I = rbs.begin();
  for (; I != rbs.end(); ++I) {
    I->getForce() += calculateSimpleGravity(I, m_g);
  };
}
