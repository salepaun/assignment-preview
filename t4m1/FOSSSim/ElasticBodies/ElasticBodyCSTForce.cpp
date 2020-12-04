#include "ElasticBodyCSTForce.h"
#include <assert.h>


#include "../Ops.h"



static inline void compGradPhi(
    Vector2s const &_Xb1,
    Vector2s const &_Xb2,
    Vector2s const &_Xb01,
    Vector2s const &_Xb02,
    Matrix2s &_GradPhi,
    Matrix2s &_AuxM)
{
  Matrix2s A, B;
  A << _Xb1.x(), _Xb2.x(), _Xb1.y(), _Xb2.y();
  B << _Xb01.x(), _Xb02.x(), _Xb01.y(), _Xb02.y();
  _AuxM = B.inverse();

  _GradPhi = A * _AuxM;

#if MY_DEBUG > 1
  D1(" 1:" << _Xb1.transpose()
      << ", 2:" << _Xb2.transpose()
      << ", 10:" << _Xb01.transpose()
      << ", 20:" << _Xb02.transpose()
      << "\nA=\n" << A
      << "\nB=\n" << B);
#endif
}


static inline void compStrain(
    Matrix2s const &_GradPhi,
    Matrix2s &_Strain)
{
  _Strain <<
    _GradPhi(0,0)*_GradPhi(0,0) + _GradPhi(1,0)*_GradPhi(1,0) - 1,
    _GradPhi(0,0)*_GradPhi(0,1) + _GradPhi(1,0)*_GradPhi(1,1),
    _GradPhi(0,0)*_GradPhi(0,1) + _GradPhi(1,0)*_GradPhi(1,1),
    _GradPhi(0,1)*_GradPhi(0,1) + _GradPhi(1,1)*_GradPhi(1,1) - 1;
}



void ElasticBodyCSTForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
  assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
  assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

  int I1, I2, I3;
  I1 = m_idx1 << 1;
  I2 = m_idx2 << 1;
  I3 = m_idx3 << 1;

  Vector2s X1 = x.segment<2>(I1);
  Vector2s X2 = x.segment<2>(I2);
  Vector2s X3 = x.segment<2>(I3);
  Vector2s Xb1 = X2-X1;
  Vector2s Xb2 = X2-X3;
  Vector2s Xb3 = X3-X1;

  Vector2s Xb01 = m_xb2 - m_xb1;
  Vector2s Xb02 = m_xb2 - m_xb3;
  Vector2s Xb03 = m_xb3 - m_xb1;

  Matrix2s GradPhi, Strain, AuxM;

  compGradPhi(Xb1, Xb3, Xb01, Xb03, GradPhi, AuxM);
  compStrain(GradPhi, Strain);

  scalar ES = 0.5 * m_youngs_modulus * areaE(Xb01, Xb03) / (1 - m_poisson_ratio*m_poisson_ratio) \
              * (Strain(0,0)*Strain(0,0) + 2*m_poisson_ratio*Strain(0,0)*Strain(1,1) + 2*(1-m_poisson_ratio)*Strain(0,1)*Strain(0,1));

#if MY_DEBUG > 0
  D1(" E=" << ES
      << ", Xb1=" << Xb1.transpose()
      << ", Xb2=" << Xb2.transpose()
      << ", Xb3=" << Xb3.transpose()
      << ", Xb01=" << Xb01.transpose()
      << ", Xb02=" << Xb02.transpose()
      << ", Xb03=" << Xb03.transpose());
#endif
}



void ElasticBodyCSTForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/2 );
  assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/2 );
  assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/2 );

  // Your code goes here!

  int I1, I2, I3;
  I1 = m_idx1 << 1;
  I2 = m_idx2 << 1;
  I3 = m_idx3 << 1;

  Vector2s X1 = x.segment<2>(I1);
  Vector2s X2 = x.segment<2>(I2);
  Vector2s X3 = x.segment<2>(I3);
  Vector2s Xb1 = X2-X1;
  Vector2s Xb2 = X2-X3;
  Vector2s Xb3 = X3-X1;

  Vector2s Xb01 = m_xb2 - m_xb1;
  Vector2s Xb02 = m_xb2 - m_xb3;
  Vector2s Xb03 = m_xb3 - m_xb1;

  Matrix2s GradPhi, GradPhi2, Strain, AuxM;
  VectorXs GradStrainXX, GradStrainYY, GradStrainXY;
  VectorXs GradA, GradB, GradC, GradD;

  compGradPhi(Xb1, Xb3, Xb01, Xb03, GradPhi, AuxM);
  compStrain(GradPhi, Strain);

  GradPhi2.Zero();
  GradPhi2 << 
    (-AuxM(0,0)-AuxM(1,0))*X1.x() + AuxM(0,0)*X2.x() + AuxM(1,0)*X3.x(),
    (-AuxM(0,1)-AuxM(1,1))*X1.x() + AuxM(0,1)*X2.x() + AuxM(1,1)*X3.x(),
    (-AuxM(0,0)-AuxM(1,0))*X1.y() + AuxM(0,0)*X2.y() + AuxM(1,0)*X3.y(),
    (-AuxM(0,1)-AuxM(1,1))*X1.y() + AuxM(0,1)*X2.y() + AuxM(1,1)*X3.y();

  GradA.resize(6);
  GradA << (-AuxM(0,0)-AuxM(1,0)), 0, AuxM(0,0), 0, AuxM(1,0), 0;

  GradB.resize(6);
  GradB << 0, (-AuxM(0,0)-AuxM(1,0)), 0, AuxM(0,0), 0, AuxM(1,0);

  GradC.resize(6);
  GradC << (-AuxM(0,1)-AuxM(1,1)), 0, AuxM(0,1), 0, AuxM(1,1), 0;

  GradD.resize(6);
  GradD << 0, (-AuxM(0,1)-AuxM(1,1)), 0, AuxM(0,1), 0, AuxM(1,1);

  scalar A,B,C,D;
  A = GradPhi(0,0); B = GradPhi(1,0); C = GradPhi(0,1); D = GradPhi(1,1);

  GradStrainXX = 2*A*GradA + 2*B*GradB;
  GradStrainYY = 2*C*GradC + 2*D*GradD;
  GradStrainXY = A*GradC + C*GradA + B*GradD + D*GradB;

  VectorXs GradE = m_youngs_modulus * areaE(Xb01, Xb03) / (1 - m_poisson_ratio*m_poisson_ratio) * ( \
                   (Strain(0,0) + m_poisson_ratio*Strain(1,1)) * GradStrainXX + \
                   (Strain(1,1) + m_poisson_ratio*Strain(0,0)) * GradStrainYY + \
                   2*(1-m_poisson_ratio) * Strain(1,0) * GradStrainXY);

  gradE.segment<2>(I1) += GradE.segment<2>(0);
  gradE.segment<2>(I2) += GradE.segment<2>(2);
  gradE.segment<2>(I3) += GradE.segment<2>(4);


#if MY_DEBUG > 0
  D1(" GradPhi=\n" << GradPhi
      << "\nGradPhi2=\n" << GradPhi2
      << "\nStrain=\n" << Strain
      << "\nGradE=" << GradE.transpose()
      << "\nGradET= 1)" << gradE.segment<2>(I1).transpose()
      << " ** 2)" << gradE.segment<2>(I2).transpose()
      << " ** 3)" << gradE.segment<2>(I3).transpose()
      << "\n ** Xb1=" << Xb1.transpose()
      << ", Xb2=" << Xb2.transpose()
      << ", Xb3=" << Xb3.transpose()
      << ", Xb01=" << Xb01.transpose()
      << ", Xb02=" << Xb02.transpose()
      << ", Xb03=" << Xb03.transpose());
#endif
}




void ElasticBodyCSTForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
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
}

void ElasticBodyCSTForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
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
