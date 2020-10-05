#include "ExplicitEuler.h"

#include <stdio.h>

#include <iostream>
#include <Eigen/Dense>

using namespace std;


static scalar g_Time = 0;


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
        const Matrix2s& _M,
        const Vector2s& _F,
        scalar _Dt)
{
    scalar DetM = _M.determinant();
    assert( DetM != 0 );
    Matrix2s N = _M.inverse();
    _V += _Dt * N * _F;
}


bool ExplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
    // Your code goes here!
    
    // Some tips on getting data from TwoDScene:
    // A vector containing all of the system's position DoFs. x0, y0, x1, y1, ...
    VectorXs& X = scene.getX();
    // A vector containing all of the system's velocity DoFs. v0, v0, v1, v1, ...
    VectorXs& V = scene.getV();
    // A vector containing the masses associated to each DoF. m0, m0, m1, m1, ...
    const VectorXs& M = scene.getM();
    // DOF
    int Size = X.size(), Num = scene.getNumParticles();
    // Kinetic Energy
    scalar Ek = 0;
    // Overall speed
    Vector2s Vt;
    // A force vector
    VectorXs F(Size);


#ifndef NDEBUG
    cout << "\nCalled:" << __FUNCTION__ \
        << ", dt=" << dt \
        << ", Size=" << Size \
        << ", Num=" << Num << endl;
#endif


    FILE *fd = NULL;
    if(g_Time) {
        fd = fopen("./kinetic_energy_output.txt", "a");
    } else {
        fd = fopen("./kinetic_energy_output.txt", "w");
        if(fd)
            fprintf(fd, "# Time\t Kinetic Energy\t Vx\t Vy\n");
        else
            perror("Failed to open kinetic energy dump");
    };


    F.setZero();
    Vt.setZero();

    // Calculate forces
    scene.accumulateGradU(F);
    
    // Zero out for fixed particles
    for(int i=0, j=0; i<Size; i+=2) {
        j = i+1;
        // Determine if the ith particle is fixed
        if(scene.isFixed(i>>1)) {
            V[i] = 0; V[j] = 0;
            F[i] = 0; F[j] = 0;
        } else {
            Vector2s Xn(X[i], X[j]);
            Vector2s Vn(V[i], V[j]);
            Matrix2s Mn; Mn << M[i], 0, 0, M[j];
            Vector2s Fn(F[i], F[j]);

            dq(Xn, Vn, dt);
            dv(Xn, Vn, Mn, Fn, dt);

            X[i] = Xn[0]; X[j] = Xn[1];
            V[i] = Vn[0]; V[j] = Vn[1];
            Vt += Vn;
        };
    };

    g_Time += dt;
    Ek = scene.computeKineticEnergy();
    if(fd) {
        fprintf(fd, "%'.2lf\t%'.5lf\t%'.5lf\t%'.5lf\n", g_Time, Ek, Vt[0], Vt[1]);
        fclose(fd);
    };

#ifndef NDEBUG
    cout \
        << "Particles:" << scene.getNumParticles() \
        << ", Size:" << X.size() \
        << ", Ek=" << Ek << endl;
    cout <<
        "\tX     :" << X.transpose() << endl;
    cout <<
        "\tV     :" << V.transpose() << endl;
    cout <<
        "\tForces:" << F.transpose() << endl;
#endif
    
    
    return true;
}
