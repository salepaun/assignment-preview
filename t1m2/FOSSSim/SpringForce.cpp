#include "SpringForce.h"

#include <Eigen/Dense>

inline static Vector2s extractVector(const VectorXs& _X, int _N)
{
  int Idx = _N << 1;
  return Vector2s(_X[Idx], _X[Idx+1]);
}

inline static void updateVector(VectorXs& _X, int _N, Vector2s const & _Xn)
{
  int Idx = _N << 1;
  _X[Idx] = _Xn.x();
  _X[Idx+1] = _Xn.y();
}

inline static scalar calculateDistance(
    Vector2s const & _Xi, \
    Vector2s const & _Xj)
{
  scalar L = 0.0;
  Vector2s Xd = _Xi - _Xj;
  // Basic
  //L = sqrt((_Xi.x-_Xj.x)*(_Xi.x-_Xj.x) + (_Xi.y-_Xj.y)*(_Xi.y-_Xj.y));
  // More general
  //L = sqrt(Xd.dot(Xd));
  // Using Eigen library
  L = Xd.norm();

  return L;
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
  E = _K/2 * (_L - _L0) * (_L - _L0);

  return E;
}

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

inline static void updateGradVector(VectorXs &_GradE, int _N, Vector2s const & _GradI)
{
  int Idx = _N << 1;
  _GradE[Idx] = _GradI.dot(Vector2s(1,0));
  _GradE[Idx+1] = _GradI.dot(Vector2s(0,1));
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
