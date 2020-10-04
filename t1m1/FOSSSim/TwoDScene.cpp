#include "TwoDScene.h"

static inline scalar kineticEnergy(const Vector2s& _V, const Matrix2s& _M)
{
    return 0.5 * (_M * _V).dot(_V);
}

scalar TwoDScene::computeKineticEnergy() const
{
  // Your code goes here!
  scalar Ek = 0;
  const VectorXs& V = getV();
  const VectorXs& M = getM();
  int Size = V.size();

  for(int i=0; i < Size; i<<=1) {
      Vector2s Vn(V[i], V[i+1]);
      Matrix2s Mn; Mn << M[i], 0, 0, M[i+1];
      Ek += kineticEnergy(Vn, Mn);
  };
  
  return Ek;
}
