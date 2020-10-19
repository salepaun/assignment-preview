#include "DragDampingForce.h"

#include "Ops.h"

#include <iostream>


using namespace std;


inline static Vector2s calculateForceGrad( \
    Vector2s const & _V, \
    scalar const &_B)
{
  return _B * _V;
}


void DragDampingForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  
  // Your code goes here!
  int Size = v.size();
  for(int i=0, n=0; i < Size; i+=2, n++) {
    Vector2s V = extractVectorIdx(x, i);
    if(!V.isZero(0)) {
      Vector2s GradEn = calculateForceGrad(V, m_b);
      updateVectorDirect(gradE, n, GradEn);
      cout << "V[" << i << ":" << n << "]=" << V.transpose() << ", E:" << GradEn.transpose() << endl;
    };
  }
}
