#include "RigidBodySpringForce.h"

#include "../Ops.h"

using namespace std;



/**
 * Calculates spring force.
 * @param _Xi
 * @param _Xj
 * @param _N a normalized vector parallel to Xij
 * @param _L0
 * @param _K
 */
static Vector2s calculateForce(
    Vector2s const &_Xi,
    Vector2s const &_Xj,
    Vector2s const &_N,
    scalar const &_L0,
    scalar const &_K)
{
  scalar L = (_Xi - _Xj).norm();
  Vector2s F = - _K * (L - _L0) * _N;

#if MY_DEBUG > 1
  D1(" F=" << F.transpose() << ", |F|=" << F.norm()
      << ", L=" << L << ", L0=" << _L0
      << ", N=" << _N.transpose());
#endif

  return F;
}


/**
 * Calculates linear and torque part of the spring force.
 * @param _F a spring force vector
 * @param _R a radius
 * @param _Rn a normal parallel to radius
 * @param _N a normal perpendicular to radius
 */
static void calculateForceTorque(
    Vector2s const &_F,
    Vector2s const &_R,
    Vector2s const &_Rn,
    Vector2s const &_Xnij,
    Vector2s const &_N,
    RigidBody &_RB)
{
  Vector2s F = _F.dot(_Xnij) * _Xnij;
  // scalar T = _F.dot(_N) * _R.norm();
  scalar T = cross2s(_R, _F);

#if MY_DEBUG > 1
  D1(" F=" << _F.transpose() << ", |F|=" << _F.norm()
      << ", Fl=" << F.transpose() << ", |Fl|=" << F.norm()
      << ", Xnij=" << _Xnij.transpose() << ", |Xnij|=" << _Xnij.norm()
      << ", R=" << _R.transpose() << ", |R|=" << _R.norm()
      << ", Torque=" << T);
#endif

  _RB.getForce() += F;
  _RB.getTorque() += T;
}



static void acquireRealWorldV(
    RigidBodySpringForce const &_Self,
    std::vector<RigidBody> const &_RBs,
    Vector2s &_Xi,
    Vector2s &_Xj,
    Vector2s &_Xij)
{
  if (_Self.firstEndpointIsFixed()) {
    _Xi = _Self.getFirstEndpoint();
  } else {
    _Xi = _Self.computeFirstEndpoint(_RBs);
  };

  if (_Self.secondEndpointIsFixed()) {
    _Xj = _Self.getSecondEndpoint();
  } else {
    _Xj = _Self.computeSecondEndpoint(_RBs);
  };

  _Xij = (_Xi - _Xj);
}



scalar RigidBodySpringForce::computePotentialEnergy( const std::vector<RigidBody>& rbs )
{
  assert( m_rb0 >= -1 ); assert( m_rb0 < (int) rbs.size() );
  assert( m_rb1 >= -1 ); assert( m_rb1 < (int) rbs.size() );
  
  // Your code goes here!
  scalar E = 0.0;

  if (isConservative())
  {
    Vector2s Xi, Xj, Xij;
    acquireRealWorldV(*this, rbs, Xi, Xj, Xij);
    scalar D = Xij.norm() - m_l0;

    E = m_k/2 * D * D;

#if MY_DEBUG > 1
    D1(", E:" << E << ", D=" << D << ", l0=" << m_l0);
#endif
  };

  return E;
}


void RigidBodySpringForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  assert( m_rb0 >= -1 ); assert( m_rb0 < (int) rbs.size() );
  assert( m_rb1 >= -1 ); assert( m_rb1 < (int) rbs.size() );
  
  // Your code goes here!
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  Vector2s Xi, Xj, Xij, Xnij;
  Vector2s Ri, Rj, N;
  int IdxRBi = getFirstRigidBody();
  int IdxRBj = getSecondRigidBody();

  acquireRealWorldV(*this, rbs, Xi, Xj, Xij);

  if (Xij.norm() != 0 || m_l0 != 0)
  {
    Xnij = Xij;
    Xnij.normalize();
    Matrix2s R90 = Matrix2s::Identity();
    R90 << 0,-1,1,0;
        
    Vector2s Force = calculateForce(Xi, Xj, Xnij, m_l0, m_k);

    if (!firstEndpointIsFixed()) {
      Ri = Xi - rbs[IdxRBi].getX();
      Vector2s Rn = Ri; Rn.normalize();
      N = R90 * Rn;
      calculateForceTorque(Force, Ri, Rn, Xnij, N, rbs[IdxRBi]);
    };
    if (!secondEndpointIsFixed()) {
      Rj = Xj - rbs[IdxRBj].getX();
      Vector2s Rn = Rj; Rn.normalize();
      N = R90 * Rn;
      calculateForceTorque(-Force, Rj, Rn, Xnij, N, rbs[IdxRBj]);
    };
  };

#if MY_DEBUG > 1
  D1(", Xi=" << Xi.transpose()
      << ", Xj=" << Xj.transpose()
      << ", Ri=" << Ri.transpose()
      << ", Rj=" << Rj.transpose())
#endif
}
