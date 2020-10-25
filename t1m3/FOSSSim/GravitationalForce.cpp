#include "GravitationalForce.h"
#include "Ops.h"

void GravitationalForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );

  // Compute the force Jacobian here!
  int i,j, Ii, Ij;
  VectorXs const &X = x;
  
  Ii = i = m_particles.first;
  Ij = j = m_particles.second;
  Ii <<= 1;
  Ij <<= 1;

  Matrix2s K = Matrix2s::Zero();
  Matrix2s Mi = m.segment<2>(Ii).asDiagonal();
  Matrix2s Mj = m.segment<2>(Ij).asDiagonal();
  Vector2s Xi = X.segment<2>(Ii);
  Vector2s Xj = X.segment<2>(Ij);
  scalar L = calculateDistance(Xi, Xj);
  Vector2s N = Xi - Xj;
  N.normalize();
  
  // Calculate 
  K = - m_G/pow(L,3) * Mi * Mj * (Matrix2s::Identity() - 3 * N * N.transpose()); 

  /* 
  Matrix4s J = Matrix4s::Zero();

  J.block<2,2>(0,0) = K;
  J.block<2,2>(2,0) = -K;
  J.block<2,2>(0,2) = -K;
  J.block<2,2>(2,2) = K;
  */

  hessE.block<2,2>(Ii,Ii) += K;
  hessE.block<2,2>(Ii,Ij) += -K;
  hessE.block<2,2>(Ij,Ii) += -K;
  hessE.block<2,2>(Ij,Ij) += K;
}

void GravitationalForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
    
  // Nothing to do.
}
