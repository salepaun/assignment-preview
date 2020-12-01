#include "ElasticBodySpringForce.h"
#include <assert.h>


#include "../Ops.h"




static inline scalar compParE(scalar const &_Alpha, scalar const &_L0)
{
  return 0.5 * _Alpha / _L0;
}


static inline scalar compParF(scalar const &_Alpha, scalar const &_L0)
{
  return _Alpha / _L0;
}



void ElasticBodySpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );

    E += compParE(m_alpha, m_l0) * pow(((x.segment<2>(m_idx1<<1) - x.segment<2>(m_idx2<<1)).norm() - m_l0), 2);
}

void ElasticBodySpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );

    // Your code goes here!
    int Idxi = m_idx1 << 1;
    int Idxj = m_idx2 << 1;
    Vector2s Xi = x.segment<2>(Idxi);
    Vector2s Xj = x.segment<2>(Idxj);
    Vector2s N = (Xi-Xj); N.normalize();

    Vector2s E = compParF(m_alpha, m_l0) * ((Xi-Xj).norm() - m_l0) * N;
    gradE.segment<2>(Idxi) += E;
    gradE.segment<2>(Idxj) += -E;
}

void ElasticBodySpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == hessE.rows() );
    assert( x.size() == hessE.cols() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    
    // Your code goes here!
    int Idxi = m_idx1 << 1;
    int Idxj = m_idx2 << 1;
    Vector2s Xi = x.segment<2>(Idxi);
    Vector2s Xj = x.segment<2>(Idxj);
    Vector2s N = (Xi-Xj); N.normalize();
    Matrix2s Id = Matrix2s::Identity();
    Matrix2s NNt = N*N.transpose();
    Matrix2s K = Matrix2s::Zero();
    scalar L = (Xi-Xj).norm();

    // Contribution from elastic component
    K = - compParF(m_alpha, m_l0) * (NNt + (L-m_l0)/L * (Id - NNt));
    updateHessianWithBlock(Idxi, Idxj, hessE, K);
}

void ElasticBodySpringForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == hessE.rows() );
    assert( x.size() == hessE.cols() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    
    // You do not need to implement this.
}
