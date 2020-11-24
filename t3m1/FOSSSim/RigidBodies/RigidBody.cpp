#include "RigidBody.h"

#include "../Ops.h"


Vector2s RigidBody::computeTotalMomentum() const
{
  // Your code goes here!
  return getM() * getV();
}

scalar RigidBody::computeCenterOfMassAngularMomentum() const
{
  // Your code goes here!DE
  return cross2s(getX(), computeTotalMomentum());
}

scalar RigidBody::computeSpinAngularMomentum() const
{
  // Your code goes here!
  return getI() * getOmega();
}

scalar RigidBody::computeCenterOfMassKineticEnergy() const
{
  // Your code goes here!
  return 0.5 * getM() * getV().dot(getV());
}

scalar RigidBody::computeSpinKineticEnergy() const
{
  // Your code goes here!
  return 0.5 * getI() * getOmega() * getOmega();
}

scalar RigidBody::computeTotalMass( const VectorXs& masses ) const
{
  // Your code goes here!
  double TotalMass = 0.0;

  for (int i=0; i < masses.size(); i++) {
    TotalMass += masses(i);
  };

#if MY_DEBUG > 0
  D1(TotalMass << " out of " << masses.size());
#endif

  return TotalMass;
}

Vector2s RigidBody::computeCenterOfMass( const VectorXs& vertices, const VectorXs& masses ) const
{
  // Your code goes here!
  double TotalMass = computeTotalMass(masses);
  Vector2s A = Vector2s::Zero();
  int Size = vertices.size();
  int VerticesNum = Size >> 1;

  for (int i=0, j=0; i < Size; i+=2, ++j) {
    A += masses(j) * vertices.segment<2>(i);
  };

  A /= TotalMass;

#if MY_DEBUG > 0
  D1(A.transpose());
#endif

  return A;
}

scalar RigidBody::computeMomentOfInertia( const VectorXs& vertices, const VectorXs& masses ) const
{
  assert( vertices.size()%2 == 0 );
  assert( 2*masses.size() == vertices.size() );
  
  // Your code goes here!
  double I = 0.0;
  int Size = vertices.size();

  for (int i=0, j=0; i < Size; i+=2, j++) {
    I += masses(j) * (vertices.segment<2>(i) - getX()).squaredNorm();
  };
  
  return I;
}
