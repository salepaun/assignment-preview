#include "ExplicitEuler.h"

#include <iostream>
#include <Eigen/Dense>

using namespace std;



static inline void dq(
        Vector2s& _Q,
        const Vector2s& _V,
        scalar _Dt)
{
    _Q += _Dt * _V;
}


static inline void dv(
        const Vector2s& _Q, 
        Vector2s& _V,
        Matrix2s& _M,
        const Vector2s& _F,
        scalar _Dt)
{
    scalar DetM = _M.determinant();
    assert( DetM != 0 );
    _M.inverse();
    _V += _Dt * _M * _F;
}


bool ExplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
    // Your code goes here!
    
    // Some tips on getting data from TwoDScene:
    // A vector containing all of the system's position DoFs. x0, y0, x1, y1, ...
    VectorXs& x = scene.getX();
    // A vector containing all of the system's velocity DoFs. v0, v0, v1, v1, ...
    VectorXs& v = scene.getV();
    // A vector containing the masses associated to each DoF. m0, m0, m1, m1, ...
    const VectorXs& m = scene.getM();
    // DOF
    int Size = x.size(), Num = scene.getNumParticles();

    // A force vector
    VectorXs f(Size);

#ifndef NDEBUG
    cout << "Called:" << __FUNCTION__ \
        << ", dt=" << dt \
        << ", Size=" << Size \
        << ", Num=" << Num << endl;
#endif

    f.setZero();

    // Calculate forces
    scene.accumulateGradU(f);
    
    // Zero out for fixed particles
    for(int i=0, j=0; i<Num; i<<=1) {
        j = i+1;
        // Determine if the ith particle is fixed
        if(scene.isFixed(i)) {
            v[i] = 0; v[j] = 0;
            f[i] = 0; f[j] = 0;
        } else {
            Vector2s Xn(x[i], x[j]);
            Vector2s Vn(v[i], v[j]);
            Matrix2s Mn; Mn << m[i], 0, 0, m[j];
            Vector2s Fn(f[i], f[j]);

            dq(Xn, Vn, dt);
            dv(Xn, Vn, Mn, Fn, dt);

            x[i] = Xn[0]; x[j] = Xn[1];
            v[i] = Vn[0]; v[j] = Vn[1];
        };
    };

#ifndef NDEBUG
    cout \
        << "Particles:" << scene.getNumParticles() \
        << ", Size:" << x.size() << endl;
    cout <<
        "\tX     :" << x << endl;
    cout <<
        "\tV     :" << v << endl;
    cout <<
        "\tForces:" << f << endl;
#endif
    
    
    return true;
}
