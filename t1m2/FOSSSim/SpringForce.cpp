#include "SpringForce.h"
#include "Ops.h"



inline static T_VecPair calculateDampingForceGrad( \
    Vector2s const & _Vi, \
    Vector2s const & _Vj, \
    Vector2s const & _Nij, \
    scalar const &_B)
{
  Vector2s GradEi, GradEj;
  GradEi = _B * _Nij.dot((_Vi - _Vj)) * _Nij;
  GradEj = -GradEi;

  return T_VecPair(GradEi, GradEj);
}

inline static T_VecPair calculateForceGrad( \
    Vector2s const & _Xi, \
    Vector2s const & _Xj, \
    Vector2s const & _Vi, \
    Vector2s const & _Vj, \
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
  T_VecPair GradDampE = calculateDampingForceGrad(_Vi, _Vj, Nij, _B);

  return T_VecPair(GradEi+GradDampE.first, GradEj+GradDampE.second);
}

inline static scalar calculateOscilatorEnergy( \
    Vector2s const & _Xi, \
    Vector2s const & _Xj, \
    scalar const &_K, \
    scalar const &_L0, \
    scalar const _L, \
    scalar const &_B)
{
  scalar E = 0.0, D = _L - _L0;
  E = _K/2 * D * D;

  return E;
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
  Vector2s Vi = extractVector(v, m_endpoints.first);
  Vector2s Vj = extractVector(v, m_endpoints.second);

  T_VecPair GradEn = calculateForceGrad(Xi, Xj, Vi, Vj, m_k, m_l0, calculateDistance(Xi, Xj), m_b);

  updateVector(gradE, m_endpoints.first, GradEn.first);
  updateVector(gradE, m_endpoints.second, GradEn.second);
}
