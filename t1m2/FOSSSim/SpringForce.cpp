#include "SpringForce.h"
#include "Ops.h"


inline static std::pair<Vector2s, Vector2s> calculateForceGrad( \
    Vector2s const & _Xi, \
    Vector2s const & _Xj, \
    scalar const &_K, \
    scalar const &_L0, \
    scalar const _L, \
    scalar const &_B)
{
  Vector2s GradEi, GradEj, Nij;
  Nij = _Xi - _Xj;
  Nij.normalize();
  GradEi = _K * (_L - _L0) * Nij;
  GradEj = -GradEi;
  return std::pair<Vector2s, Vector2s>(GradEi, GradEj);
}


void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
  Vector2s Xi = extractVector(x, m_endpoints.first);
  Vector2s Xj = extractVector(x, m_endpoints.second);

  E += calculateOscilatorEnergy(Xi, Xj, m_k, m_l0, calculateDistance(Xi, Xj), m_b);
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
  Vector2s Xi = extractVector(x, m_endpoints.first);
  Vector2s Xj = extractVector(x, m_endpoints.second);

  std::pair<Vector2s, Vector2s> GradEi = 
    calculateForceGrad(Xi, Xj, m_k, m_l0, calculateDistance(Xi, Xj), m_b);

  updateGradVector(gradE, m_endpoints.first, GradEi.first);
  updateGradVector(gradE, m_endpoints.second, GradEi.second);
}
