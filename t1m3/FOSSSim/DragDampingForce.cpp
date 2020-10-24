#include "DragDampingForce.h"

#include <iostream>


using namespace std;



void DragDampingForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );

  // Nothing to do.
}

void DragDampingForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  
  // Compute the force Jacobian here!
  int Size = x.size();
  for(int i=0; i < Size; i++) {
    hessE(i,i) -= m_b;
  };

  cout << "Called:" << __FUNCTION__ << ": H:" << hessE << endl;
}
