#include "GravitationalForce.h"
#include "Ops.h"

#include <iostream>

using namespace std;

inline static scalar calculateGravityEnergy( \
    Vector2s const & _Xi, \
    Vector2s const & _Xj, \
    Vector2s const & _Mi, \
    Vector2s const & _Mj, \
    scalar const &_G, \
    scalar const _L)
{
  scalar E = -_G/_L * _Mi.dot(_Mj);
  return E;
}

inline static std::pair<Vector2s, Vector2s> calculateForceGrad( \
    Vector2s const & _Xi, \
    Vector2s const & _Xj, \
    Matrix2s const & _Mi, \
    Matrix2s const & _Mj, \
    scalar const &_G, \
    scalar const _L2)
{
  Vector2s GradEi, GradEj, Nij;
  Nij = _Xi - _Xj;
  Nij.normalize();
  GradEi = _G/_L2 * _Mi * _Mj * Nij;
  GradEj = -GradEi;
  return std::pair<Vector2s, Vector2s>(GradEi, GradEj);
}


void GravitationalForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
  Vector2s Xi = extractVector(x, m_particles.first);
  Vector2s Xj = extractVector(x, m_particles.second);
  Vector2s Mi = extractVector(m, m_particles.first);
  Vector2s Mj = extractVector(m, m_particles.second);
  E += calculateGravityEnergy(Xi, Xj, Mi, Mj, m_G, calculateDistance(Xi, Xj));

#ifndef NDEBUG
  /* std::cout << "Mi:" << Mi.transpose() << "\n" \
    << "Mj:" << Mj.transpose() << std::endl; */
#endif
}

void GravitationalForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
  Vector2s Xi = extractVector(x, m_particles.first);
  Vector2s Xj = extractVector(x, m_particles.second);
  Matrix2s Mi = extractMassMx(m, m_particles.first);
  Matrix2s Mj = extractMassMx(m, m_particles.second);
  std::pair<Vector2s, Vector2s> GradEn = calculateForceGrad(Xi, Xj, Mi, Mj, m_G, calculateSquaredDistance(Xi, Xj));

  updateVector(gradE, m_particles.first, GradEn.first);
  updateVector(gradE, m_particles.second, GradEn.second);

#ifndef NDEBUG
  /*
  Vector2s Miv = extractVector(m, m_particles.first);
  Vector2s Mjv = extractVector(m, m_particles.second);
  std::cout << "Mi:" << Miv.transpose() << "\n" \
    << "Mj:" << Mjv.transpose() << std::endl;
  */
#endif
}
