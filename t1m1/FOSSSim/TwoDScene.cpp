#include "TwoDScene.h"

#include <stdio.h>

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

    for(int i=0,j=0; i<Size; i+=2) {
        j=i+1;
        Vector2s Vn(V[i], V[j]);
        Matrix2s Mn; Mn << M[i], 0, 0, M[j];
        Ek += kineticEnergy(Vn, Mn);
    };

    return Ek;
}
