#include "ElasticBodyBendingForce.h"
#include <assert.h>


#include "../Ops.h"



static inline scalar compParE(scalar const &_Alpha)
{
    return 0.5 * _Alpha;
}


static inline scalar compParF(scalar const &_Alpha, scalar const &_Theta0)
{
    return 0.5 * _Alpha;
}

static inline scalar compTheta(
        Vector2s const &_Eij,
        Vector2s const &_Ejk)
{
    return atan2(cross2s(_Eij, _Ejk), _Eij.dot(_Ejk));
}


static inline scalar compGradEPar(
        scalar const &_Alpha,
        scalar const &_Theta0,
        scalar const &_Eij0,
        scalar const &_Ejk0,
        Vector2s const &_Eij,
        Vector2s const &_Ejk)
{
    scalar Theta = compTheta(_Eij, _Ejk);
    scalar A = _Alpha / (_Eij0 + _Ejk0) * (Theta - _Theta0);

#if MY_DEBUG > 1
    D1("A=" << A
            << ", Alpha=" << _Alpha
            << ", Theta0=" << _Theta0
            << ", Theta=" << Theta);
#endif

    return A;
}


static inline void compGradAtan2(
        Vector2s const &_Eij,
        Vector2s const &_Ejk,
        Vector2s &_GradThetaPar)
{
    scalar X, Y, B;
    Y = cross2s(_Eij, _Ejk);
    X = _Eij.dot(_Ejk);
    B = X*X + Y*Y;
    _GradThetaPar << (-Y / B), (X / B);
    //_GradThetaPar << (X / B), (-Y / B);

#if MY_DEBUG > 1
    D1("Eij=" << _Eij.transpose()
            << ", Ejk=" << _Ejk.transpose()
            << ", Y=" << Y
            << ", X=" << X
            << ", B=" << B);
#endif
}



static inline void compGradX(
        Vector2s const &_Eij,
        Vector2s const &_Ejk,
        Matrix2s const &_R,
        MatrixXs &_GradX)
{
    _GradX.resize(3,2);
    _GradX.row(0) = _Ejk;
    _GradX.row(1) = _Eij - _Ejk;
    _GradX.row(2) = -_Eij;
}



static inline void compGradY(
        Vector2s const &_Eij,
        Vector2s const &_Ejk,
        Matrix2s const &_R,
        MatrixXs &_GradY)
{
    _GradY.resize(3,2);
    _GradY.row(0) = _R * _Ejk;
    _GradY.row(1) = _R * (-_Eij - _Ejk);
    _GradY.row(2) = _R * _Eij;
}



void ElasticBodyBendingForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

    int Ii, Ij, Ik;
    Vector2s Xi, Xj, Xk, Eij, Ejk;

    Ii = m_idx1 << 1; Ij = m_idx2 << 1; Ik = m_idx3 << 1;
    Xi = x.segment<2>(Ii);
    Xj = x.segment<2>(Ij);
    Xk = x.segment<2>(Ik);
    Eij = Xi - Xj;
    Ejk = Xj - Xk;

    scalar EE = compParE(m_alpha) * 1 / (Eij.norm() + Ejk.norm()) * pow((compTheta(Eij, Ejk) - m_theta0), 2);

#if MY_DEBUG > 0
    D1(" Elastic Energy=" << EE);
#endif
    E += EE;
}



void ElasticBodyBendingForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

    // Your code goes here!

    int Ii, Ij, Ik;
    Vector2s Xi, Xj, Xk, Eij, Ejk, Eik, Nij, Njk, Nik;
    Matrix2s I = Matrix2s::Identity();
    Matrix2s R90; R90 << 0,-1,1,0;
    Matrix2s R270; R270 << 0,1,-1,0;
    Matrix2s Id = Matrix2s::Identity();
    Vector2s GradEPar, GradEiPar, GradEjPar, GradEkPar, GradEi, GradEj, GradEk;

    Ii = m_idx1 << 1; Ij = m_idx2 << 1; Ik = m_idx3 << 1;
    Xi = x.segment<2>(Ii);
    Xj = x.segment<2>(Ij);
    Xk = x.segment<2>(Ik);
    Eij = Xi - Xj;
    Ejk = Xj - Xk;
    Eik = Xi - Xk;
    Nij = Eij;
    Njk = Ejk;
    Nik = Eik;
    Nij.normalize();
    Njk.normalize();
    Nik.normalize();

    scalar A = compGradEPar(m_alpha, m_theta0, m_eb1n, m_eb2n, Eij, Ejk);

    compGradAtan2(Eij, Ejk, GradEPar);

    /*
    GradEiPar = GradEPar.x() * Ejk + GradEPar.y() * R90 * Ejk;
    GradEi = A * GradEiPar;

    GradEjPar = GradEPar.x() * (Eij - Ejk) + GradEPar.y() * R90 *(-Ejk - Eij);
    GradEj = A * GradEjPar;

    GradEkPar = GradEPar.x() * (-Eij) + GradEPar.y() * R90 * (Eij);
    GradEk = A * GradEkPar;
    */

    /*
    GradEiPar = Vector2s(GradEPar.x() * Nij.dot(Ejk), GradEPar.y() * cross2s(Nij, Ejk));
    GradEi = A * GradEiPar;

    GradEjPar = Vector2s(GradEPar.x() * (Njk.dot(Eij) - Nij.dot(Ejk)), GradEPar.y() * (-cross2s(Nij,Ejk) - cross2s(Njk, Eij)));
    GradEj = A * GradEjPar;

    GradEkPar = Vector2s(GradEPar.x() * -Njk.dot(Eij), GradEPar.y() * -cross2s(Njk, Eij));
    GradEk = A * GradEkPar;
    */

    /*
    GradEiPar = GradEPar.x() * Ejk + GradEPar.y() * R270 * Ejk;
    GradEi = A * GradEiPar;

    GradEjPar = GradEPar.x() * (Eij - Ejk) + GradEPar.y() * R270 *(-Ejk - Eij);
    GradEj = A * GradEjPar;

    GradEkPar = GradEPar.x() * (-Eij) + GradEPar.y() * R270 * (Eij);
    GradEk = A * GradEkPar;
    */

    MatrixXs GradE(3,2), GradQX(3,2), GradQY(3,2), GradQ(3,2);

    compGradX(Eij, Ejk, R270, GradQX);
    compGradY(Eij, Ejk, R270, GradQY);
    GradQ = GradEPar.x() * GradQX + GradEPar.y() * GradQY;
    GradE = A * GradQ;

    gradE.segment<2>(Ii) += GradE.row(0);
    gradE.segment<2>(Ij) += GradE.row(1);
    gradE.segment<2>(Ik) += GradE.row(2);

#if MY_DEBUG > 0
    D1(" Par=" << A
            << "\n ** Eij=" << Eij.transpose() << ", Ejk=" << Ejk.transpose() << ", Eik=" << Eik.transpose()
            << "\n ** GradEPar = " << GradEiPar.transpose()
            << "\n ** GradQ:\n" << GradQ
            << "\n ** GradE:\n" << GradE
            << "\n ** GradEiT= " << gradE.segment<2>(Ii).transpose()
            << " ** GradEjT= " << gradE.segment<2>(Ij).transpose()
            << " ** GradEkT= " << gradE.segment<2>(Ik).transpose());
#endif 
}




void ElasticBodyBendingForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == hessE.rows() );
    assert( x.size() == hessE.cols() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

    // Your code goes here!

    int Ii, Ij, Ik;
    Vector2s Xi, Xj, Xk, Eij, Ejk;

    Ii = m_idx1 << 1; Ij = m_idx2 << 1; Ik = m_idx3 << 1;
    Xi = x.segment<2>(Ii);
    Xj = x.segment<2>(Ij);
    Xk = x.segment<2>(Ik);
    Eij = Xi - Xj;
    Ejk = Xj - Xk;

}

void ElasticBodyBendingForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == hessE.rows() );
    assert( x.size() == hessE.cols() );
    assert( x.size()%2 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

    // You do not need to implement this.
}
