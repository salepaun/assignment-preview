#include "SimpleGravityForce.h"
#include <iostream>

using namespace std;

void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%2 == 0 );

    // Your code goes here!
    int Size=x.size(), Num=Size>>1;

#ifndef NDEBUG
    cout << "Called:" << __FUNCTION__ \
        << ", Size=" << Size \
        << ", Num=" << Num \
        << ", E=" << E << endl;
#endif
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );

    // Your code goes here!
    int Size=x.size(), Num=Size>>1;

#ifndef NDEBUG
    cout << "Called:" << __FUNCTION__ \
        << ", Size=" << Size \
        << ", Num=" << Num \
        << ", gradE.size=" << gradE.size() << endl;
#endif

    scalar LocGradE = 0;

    for(int i=0,j=0; i<Size; i+=2) {
        j=i+1;
        Matrix2s Mn; Mn << m[i], 0, 0, m[j];
        Vector2s& G = m_gravity;
        Vector2s F = Mn * G;
        gradE[i] = F[0]; gradE[j] = F[1];
    };
}
