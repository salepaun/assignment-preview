#include "RigidBodyExplicitEuler.h"

#include "../Ops.h"

bool RigidBodyExplicitEuler::stepScene( RigidBodyScene& scene, scalar dt )
{
  // Clear any previously computed forces in the rigid bodies
  std::vector<RigidBody>& rbs = scene.getRigidBodies();
  for( std::vector<RigidBody>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].getForce().setZero();
  for( std::vector<RigidBody>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].getTorque() = 0.0;
  
  // Add each force's contribution to each rigid body using previous step's state
  std::vector<RigidBodyForce*>& frcs = scene.getForces();
  for( std::vector<RigidBodyForce*>::size_type i = 0; i < frcs.size(); ++i ) frcs[i]->computeForceAndTorque(rbs);
  
  // Useful method:
  //   RigidBody::getX()
  //   RigidBody::getV()
  //   RigidBody::getTheta()
  //   RigidBody::getOmega()
  //   RigidBody::getForce()
  //   RigidBody::getTorque()
  
  // For each rigid body
  T_RBodyCntr::iterator I = rbs.begin();
  for(; I != rbs.end(); ++I)
  {
    // Your code goes here!
    Vector2s V = I->getV();
    scalar Omega = I->getOmega();

    I->getV() += dt * I->getForce() / I->getM();
    I->getOmega() += dt * I->getTorque() / I->getI();

    I->getX() += dt * V;
    I->getTheta() += dt * Omega;
  };
  
  for( std::vector<RigidBodyForce*>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].updateDerivedQuantities();
  
  return true;
}
