#include "SpringForce.h"
#include "Ops.h"

void SpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!
  int i,j, Ii, Ij;
  VectorXs const &X = x;
  VectorXs const &V = v;
  
  Ii = i = m_endpoints.first;
  Ij = j = m_endpoints.second;
  Ii <<= 1;
  Ij <<= 1;

  Matrix2s K = Matrix2s::Zero();
  Vector2s Xi = X.segment<2>(Ii);
  Vector2s Xj = X.segment<2>(Ij);
  Vector2s Vi = V.segment<2>(Ii);
  Vector2s Vj = V.segment<2>(Ij);
  Vector2s Vij = Vi - Vj;
  scalar L = calculateDistance(Xi, Xj);
  Vector2s N = Xi - Xj;
  Matrix2s Id = Matrix2s::Identity();
  N.normalize();
  Matrix2s NNt = N*N.transpose();
  
  // Contribution from elastic component
  K = -m_k*(NNt + (L-m_l0)/L * (Id - NNt));

  // Contribution from damping
  K += -m_b/L * (N.dot(Vij)*Id + N*Vij.transpose()) * (Id - NNt);

  updateHessianWithJacobianBlock(Ii, Ij, hessE, K);
}

void SpringForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!
  int i,j, Ii, Ij;
  VectorXs const &X = x;
  
  Ii = i = m_endpoints.first;
  Ij = j = m_endpoints.second;
  Ii <<= 1;
  Ij <<= 1;

  Matrix2s K = Matrix2s::Zero();
  Vector2s Xi = X.segment<2>(Ii);
  Vector2s Xj = X.segment<2>(Ij);
  Vector2s N = Xi - Xj;
  N.normalize();
  Matrix2s NNt = N*N.transpose();

  // Contribution from damping
  K = -m_b * NNt;

  updateHessianWithJacobianBlock(Ii, Ij, hessE, K);
}
