#include "SimpleCollisionHandler.h"
#include <iostream>
#include <set>

#include "MathDefs.h"

#include <vector>
// #include <lapacke.h>
#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <iomanip>

#include "Ops.h"



using namespace std;


#define IS_INF(a) (isinf(a(0,0)) )
#define IS_NOT_INF(a) (!isinf(a(0,0)))
#define CHECK_V(a) (!isinf(a(0,0)) && !isnan(a(0,0)))
#define ASSERT_MASS_CHECK(a) assert(a(0,0) == a(1,1))



static void calculateVertexVertexResponse(
    int _I1,
    int _I2,
    Vector2s const &_V1,
    Vector2s const &_V2,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    bool _Fixed1,
    bool _Fixed2,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);

    _V1p = _V2p = Vector2s::Zero();

    if (_Fixed1 || _Fixed2) {
        if (!_Fixed1) {
            _V1p = _V1 + _Vaux * _N;
            _V.segment<2>(_I1) = _V1p;
        } else if (!_Fixed2) {
            _V2p = _V2 - _Vaux * _N;
            _V.segment<2>(_I2) = _V2p;
        };
    }
    else
    {
        Matrix2s M12 = _M2 + _M1;
        Matrix2s M12i = M12.inverse();
        _V1p = _V1 + _M2 * M12i * _Vaux * _N;
        _V2p = _V2 - _M1 * M12i * _Vaux * _N;

        _V.segment<2>(_I1) = _V1p;
        _V.segment<2>(_I2) = _V2p;
    };
}



static void calculateVertexVertexScalarMResponse(
    int _I1,
    int _I2,
    Vector2s const &_V1,
    Vector2s const &_V2,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    bool _Fixed1,
    bool _Fixed2,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);
    scalar M1, M2;
    M1 = _M1(0,0);
    M2 = _M2(0,0);

    _V1p = _V2p = Vector2s::Zero();

    if (_Fixed1 || _Fixed2) {
        if (!_Fixed1) {
            _V1p = _V1 + _Vaux * _N;
            _V.segment<2>(_I1) = _V1p;
        } else if (!_Fixed2) {
            _V2p = _V2 - _Vaux * _N;
            _V.segment<2>(_I2) = _V2p;
        };
    }
    else
    {
        scalar M12 = M2 + M1;
        scalar M12i = 1/M12;
        _V1p = _V1 + M2 * M12i * _Vaux * _N;
        _V2p = _V2 - M1 * M12i * _Vaux * _N;

        _V.segment<2>(_I1) = _V1p;
        _V.segment<2>(_I2) = _V2p;
    };
}



static void calculateVertexVertexScalarMInfResponse(
    int _I1,
    int _I2,
    Vector2s const &_V1,
    Vector2s const &_V2,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    bool _Fixed1,
    bool _Fixed2,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);
    scalar M1, M2;
    M1 = _Fixed1 ? numeric_limits<scalar>::infinity() : _M1(0,0);
    M2 = _Fixed2 ? numeric_limits<scalar>::infinity() : _M2(0,0);

    _V1p = _V2p = Vector2s::Zero();

    _V1p = _V1 + (_Vaux / (1 + M1/M2)) * _N;
    _V2p = _V2 - (_Vaux / (M2/M1 + 1)) * _N;

    if (CHECK_V(_V1p)) _V.segment<2>(_I1) = _V1p;
    if (CHECK_V(_V2p)) _V.segment<2>(_I2) = _V2p;
}



static void calculateVertexEdgeResponse(
    int _I1,
    int _I2,
    int _I3,
    Vector2s const &_V1,
    Vector2s const &_V2,
    Vector2s const &_V3,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    Matrix2s const &_M3,
    bool _Fixed1,
    bool _Fixed2,
    bool _Fixed3,
    scalar const &_sAlpha,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p,
    Vector2s &_V3p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);
    ASSERT_MASS_CHECK(_M3);

    _V1p = _V2p = _V3p = Vector2s::Zero();

    scalar sBeta = 1 - _sAlpha;
    scalar sAlpha2 = _sAlpha*_sAlpha;
    scalar sBeta2 = sBeta*sBeta;

    if (_Fixed1 || _Fixed2 || _Fixed3) {
        if (!_Fixed1 && !_Fixed2)
        {
            Matrix2s M12 = _M2 + _M1*sBeta2;
            Matrix2s M12i = M12.inverse();
            _V1p = _V1 + (_M2 * M12i * _Vaux) * _N;
            _V2p = _V2 - (_M1 * M12i * _Vaux * sBeta) * _N;

            _V.segment<2>(_I1) = _V1p;
            _V.segment<2>(_I2) = _V2p;
        }
        else if (!_Fixed1 && !_Fixed3)
        {
            Matrix2s M13 = _M3 + _M1*sAlpha2;
            Matrix2s M13i = M13.inverse();
            _V1p = _V1 + (_M3 * M13i * _Vaux) * _N;
            _V3p = _V3 - (_M1 * M13i * _Vaux * _sAlpha) * _N;

            _V.segment<2>(_I1) = _V1p;
            _V.segment<2>(_I3) = _V3p;

            //_V1p = _V1 + _Vaux / (1 + sBeta2*m1/m2 + sAlpha2*m1/m3);
        }
        else if (!_Fixed2 && !_Fixed3)
        {
            Matrix2s M23 = _M3*sBeta2 + _M2*sAlpha2;
            Matrix2s M23i = M23.inverse();
            _V2p = _V2 - (_M3 * M23i * _Vaux * sBeta) * _N;
            _V3p = _V3 - (_M2 * M23i * _Vaux * _sAlpha) * _N;

            _V.segment<2>(_I2) = _V2p;
            _V.segment<2>(_I3) = _V3p;
        }
        else if (!_Fixed1)
        {
            _V1p = _V1 + _Vaux * _N;
            _V.segment<2>(_I1) = _V1p;
        }
        else if (!_Fixed2)
        {
            _V2p = _V2 - (_Vaux / sBeta) * _N;
            _V.segment<2>(_I2) = _V2p;
        }
        else if (!_Fixed3)
        {
            _V3p = _V3 - (_Vaux / _sAlpha) * _N;
            _V.segment<2>(_I3) = _V3p;
        };
    }
    else
    {
        // For mass matrixes with only diagnoal _M1 * _M2 = _M2 * _M1
        Matrix2s M12 = _M1 * _M2;
        Matrix2s M13 = _M1 * _M3;
        Matrix2s M23 = _M2 * _M3;
        Matrix2s Maux = M23 + _M1 * (sBeta2*_M3 + sAlpha2*_M2);
        Matrix2s MauxI = Maux.inverse();

        _V1p = _V1 + (M23 * MauxI * _Vaux) * _N;
        _V2p = _V2 - (M13 * MauxI * _Vaux * sBeta) * _N;
        _V3p = _V3 - (M12 * MauxI * _Vaux * _sAlpha) * _N;

        _V.segment<2>(_I1) = _V1p;
        _V.segment<2>(_I2) = _V2p;
        _V.segment<2>(_I3) = _V3p;
    };
}



static void calculateVertexEdgeScalarMResponse(
    int _I1,
    int _I2,
    int _I3,
    Vector2s const &_V1,
    Vector2s const &_V2,
    Vector2s const &_V3,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    Matrix2s const &_M3,
    bool _Fixed1,
    bool _Fixed2,
    bool _Fixed3,
    scalar const &_sAlpha,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p,
    Vector2s &_V3p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);
    ASSERT_MASS_CHECK(_M3);

    _V1p = _V2p = _V3p = Vector2s::Zero();

    scalar sBeta = 1 - _sAlpha;
    scalar sAlpha2 = _sAlpha*_sAlpha;
    scalar sBeta2 = sBeta*sBeta;
    scalar M1, M2, M3;
    M1 = _M1(0,0);
    M2 = _M2(0,0);
    M3 = _M3(0,0);

    if (_Fixed1 || _Fixed2 || _Fixed3) {
        if (!_Fixed1 && !_Fixed2)
        {
            scalar M12 = M2 + M1*sBeta2;
            scalar M12i = 1/M12;
            _V1p = _V1 + (M2 * M12i * _Vaux) * _N;
            _V2p = _V2 - (M1 * M12i * _Vaux * sBeta) * _N;

            _V.segment<2>(_I1) = _V1p;
            _V.segment<2>(_I2) = _V2p;
        }
        else if (!_Fixed1 && !_Fixed3)
        {
            scalar M13 = M3 + M1*sAlpha2;
            scalar M13i = 1/M13;
            _V1p = _V1 + (M3 * M13i * _Vaux) * _N;
            _V3p = _V3 - (M1 * M13i * _Vaux * _sAlpha) * _N;

            _V.segment<2>(_I1) = _V1p;
            _V.segment<2>(_I3) = _V3p;

            //_V1p = _V1 + _Vaux / (1 + sBeta2*m1/m2 + sAlpha2*m1/m3);
        }
        else if (!_Fixed2 && !_Fixed3)
        {
            scalar M23 = M3*sBeta2 + M2*sAlpha2;
            scalar M23i = 1/M23;
            _V2p = _V2 - (M3 * M23i * _Vaux * sBeta) * _N;
            _V3p = _V3 - (M2 * M23i * _Vaux * _sAlpha) * _N;

            _V.segment<2>(_I2) = _V2p;
            _V.segment<2>(_I3) = _V3p;
        }
        else if (!_Fixed1)
        {
            _V1p = _V1 + (_Vaux) * _N;
            _V.segment<2>(_I1) = _V1p;
        }
        else if (!_Fixed2)
        {
            _V2p = _V2 - (_Vaux / sBeta) * _N;
            _V.segment<2>(_I2) = _V2p;
        }
        else if (!_Fixed3)
        {
            _V3p = _V3 - (_Vaux / _sAlpha) * _N;
            _V.segment<2>(_I3) = _V3p;
        };
    }
    else
    {
        // For mass matrixes with only diagnoal M1 * M2 = M2 * M1
        scalar M12 = M1 * M2;
        scalar M13 = M1 * M3;
        scalar M23 = M2 * M3;
        scalar Maux = M23 + M1 * (sBeta2*M3 + sAlpha2*M2);
        scalar MauxI = 1/Maux;

        _V1p = _V1 + (M23 * MauxI * _Vaux) * _N;
        _V2p = _V2 - (M13 * MauxI * _Vaux * sBeta) * _N;
        _V3p = _V3 - (M12 * MauxI * _Vaux * _sAlpha) * _N;

        _V.segment<2>(_I1) = _V1p;
        _V.segment<2>(_I2) = _V2p;
        _V.segment<2>(_I3) = _V3p;
    };
}



static void calculateVertexEdgeScalarMInfResponse(
    int _I1,
    int _I2,
    int _I3,
    Vector2s const &_V1,
    Vector2s const &_V2,
    Vector2s const &_V3,
    scalar const &_Vaux,
    Vector2s const &_N,
    Matrix2s const &_M1,
    Matrix2s const &_M2,
    Matrix2s const &_M3,
    bool _Fixed1,
    bool _Fixed2,
    bool _Fixed3,
    scalar const &_sAlpha,
    VectorXs &_V,
    Vector2s &_V1p,
    Vector2s &_V2p,
    Vector2s &_V3p
    )
{
    ASSERT_MASS_CHECK(_M1);
    ASSERT_MASS_CHECK(_M2);
    ASSERT_MASS_CHECK(_M3);

    _V1p = _V2p = _V3p = Vector2s::Zero();

    scalar sBeta = 1 - _sAlpha;
    scalar sAlpha2 = _sAlpha*_sAlpha;
    scalar sBeta2 = sBeta*sBeta;
    scalar M1, M2, M3;
    M1 = _Fixed1 ? numeric_limits<scalar>::infinity() : _M1(0,0);
    M2 = _Fixed2 ? numeric_limits<scalar>::infinity() : _M2(0,0);
    M3 = _Fixed3 ? numeric_limits<scalar>::infinity() : _M3(0,0);

    _V1p = _V1 + (_Vaux / (1 + sBeta2*M1/M2 + sAlpha2*M1/M3)) * _N;
    _V2p = _V2 - (_Vaux*sBeta / (M2/M1 + sBeta2 + sAlpha2*M2/M3)) * _N;
    _V3p = _V3 - (_Vaux*_sAlpha / (M3/M1 + sBeta2*M3/M2 + sAlpha2)) * _N;

    if (CHECK_V(_V1p)) _V.segment<2>(_I1) = _V1p;
    if (CHECK_V(_V2p)) _V.segment<2>(_I2) = _V2p;
    if (CHECK_V(_V3p)) _V.segment<2>(_I3) = _V3p;
}

    

// BEGIN STUDENT CODE //


// Detects whether two particles are overlapping (including the radii of each)
// and approaching.
// If the two particles overlap and are approaching, returns true and sets 
// the vector n to be the vector between the first and second particle.
// Inputs:
//   scene: The scene data structure. The positions and radii of the particles
//          can be obtained from here.
//   idx1:  The index of the first particle. (Ie, the degrees of freedom
//          corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//          scene.getX().
//   idx2:  The index of the second particle.
// Outputs:
//   n: The vector between the two particles.
//   Returns true if the two particles overlap and are approaching.
bool SimpleCollisionHandler::detectParticleParticle(TwoDScene &scene, int idx1, int idx2, Vector2s &n)
{
    // Your code goes here!
    int I1 = idx1 << 1;
    int I2 = idx2 << 1;
    bool bTouch, bApproach;

    VectorXs const &X = scene.getX();
    VectorXs const &V = scene.getV();

    vector<scalar> const & Radii = scene.getRadii();

    Vector2s X1 = X.segment<2>(I1);
    Vector2s X2 = X.segment<2>(I2);

    Vector2s V1 = V.segment<2>(I1);
    Vector2s V2 = V.segment<2>(I2);

    scalar R1 = Radii[idx1];
    scalar R2 = Radii[idx2];

    Vector2s X1X2AxN = X1 - X2;
    // output according to description (not normalized!)
    n = -X1X2AxN;
    scalar D = X1X2AxN.norm();
    X1X2AxN.normalize();

    bTouch = D < R1 + R2;

    if(bTouch) {
        Vector2s Vd = V1 - V2;
        Vector2s N = n;
        N.normalize();
        bApproach = Vd.dot(N) > 0;
    } else {
        bApproach = false;
    };

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(10) \
        << "\n  X1:" << X1.transpose() << ", X2:" << X2.transpose() \
        << ", D:" << D << ", R1:" << R1 << ", R2:" << R2 \
        << ": T:" << bTouch << ", A:" << bApproach << ", N=(" << n.transpose() << ")" \
        << endl;
#endif

    return bTouch && bApproach;
}

// Detects whether a particle and an edge are overlapping (including the radii 
// of both) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the edge.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge. (Ie, the indices of particle with index e are
//          scene.getEdges()[e].first and scene.getEdges()[e].second.)
// Outputs:
//   n: The shortest vector between the particle and the edge.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleEdge(TwoDScene &scene, int vidx, int eidx, Vector2s &n)
{

    // Your code goes here!
    bool bTouch, bApproach;
    pair<int, int> const &E = scene.getEdges()[eidx];

    int I1 = vidx << 1;
    int I2 = E.first << 1;
    int I3 = E.second << 1;

    VectorXs const &X = scene.getX();
    VectorXs const &V = scene.getV();

    vector<scalar> const & Radii = scene.getRadii();

    Vector2s X1 = X.segment<2>(I1);
    Vector2s X2 = X.segment<2>(I2);
    Vector2s X3 = X.segment<2>(I3);

    Vector2s V1 = V.segment<2>(I1);
    Vector2s V2 = V.segment<2>(I2);
    Vector2s V3 = V.segment<2>(I3);

    scalar R1 = Radii[vidx];
    scalar R2 = Radii[E.first];
    scalar R3 = Radii[E.second];
    scalar Re = scene.getEdgeRadii()[eidx];
    scalar R = Re; // seems like oracle uses edge, max(max(R2,R3),Re); //max3(R2,R3,Re);


    // alpha parameter (edge drop point absolute distance)
    Vector2s Vx;
    scalar sAlpha;

    findAlpha(X1, X2, X3, V1, V2, V3, n, Vx, sAlpha);

    scalar D = n.norm();

    bTouch = D < R1 + R;

    if(bTouch) {
        Vector2s Vd = V1 - Vx;
        Vector2s N = n;
        N.normalize();
        bApproach = Vd.dot(N) > 0;
    } else {
        bApproach = false;
    };

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(10) \
        << "\n  X1:" << X1.transpose() << ", X2:" << X2.transpose() << ", X3:" << X3.transpose() \
        << ", D:" << D << ", R1:" << R1 << ", R2:" << R2 << ", R3:" << R3 << ", Re:" << Re << ", R:" << R \
        << ", Alpha:" << sAlpha << ", Xx:" << n.transpose() \
        << ", Vx:" << Vx.transpose() << ", V1:" << V1.transpose() \
        << ": T:" << bTouch << ", A:" << bApproach << ", N=(" << n.transpose() << ")" \
        << endl;
#endif

    return bTouch && bApproach;
}

// Detects whether a particle and a half-plane are overlapping (including the 
// radius of the particle) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the half-plane.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the halfplane. The vectors (px, py) and (nx, ny) can
//          be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n: The shortest vector between the particle and the half-plane.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleHalfplane(TwoDScene &scene, int vidx, int pidx, Vector2s &n)
{
    // Your code goes here!
    int I1 = vidx << 1;
    bool bTouch, bApproach;

    VectorXs const &X = scene.getX();
    VectorXs const &V = scene.getV();

    VectorXs Xp = scene.getHalfplane(pidx).first;
    VectorXs Np = scene.getHalfplane(pidx).second;

    Vector2s X1 = X.segment<2>(I1);
    Vector2s V1 = V.segment<2>(I1);
    scalar R1 = scene.getRadii()[vidx];

    Np.normalize(); // just in case
    Vector2s X1PxAx = X1 - Xp;
    Vector2s X1PAx = X1PxAx.dot(Np) * Np;

    // output according to description (not normalized!)
    n = -X1PAx;

    scalar D = X1PAx.norm();

    bTouch = D < R1;

    if(bTouch) {
        Vector2s N = n;
        N.normalize();
        bApproach = V1.dot(N) > 0;
    } else {
        bApproach = false;
    };

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(10) \
        << ", X1:" << X1.transpose() << ", Xp:" << Xp.transpose() \
        << ", D:" << D << ", R1:" << R1 \
        << ": T:" << bTouch << ", A:" << bApproach << ", N=(" << n.transpose() << ")" \
        << endl;
#endif

    return bTouch && bApproach;
}


// Responds to a collision detected between two particles by applying an impulse
// to the velocities of each one.
// You can get the COR of the simulation by calling getCOR().
// Inputs:
//   scene: The scene data structure.
//   idx1:  The index of the first particle.
//   idx2:  The index of the second particle.
//   n:     The vector between the first and second particle.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleParticle(TwoDScene &scene, int idx1, int idx2, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    VectorXs &V = scene.getV();

    // Your code goes here!
    int I1 = idx1 << 1;
    int I2 = idx2 << 1;
    Vector2s V1 = V.segment<2>(I1);
    Vector2s V2 = V.segment<2>(I2);
    Vector2s V1p, V2p;
    Vector2s N = n;

    N.normalize();
    scalar C = (1.0 + getCOR());
    scalar Vaux = C * (V2-V1).dot(N);

    /*
       Vector2s Mv1 = M.segment<2>(I1);
       Vector2s Mv2 = M.segment<2>(I2);
       Matrix2s M12 = Vector2s(Mv1.x() / Mv2.x(), Mv1.y() / Mv2.y()).asDiagonal();
       Matrix2s M21 = Vector2s(Mv2.x() / Mv1.x(), Mv2.y() / Mv1.y()).asDiagonal();
       M12 += Matrix2s::Identity(); M12.inverse();
       M21 += Matrix2s::Identity(); M21.inverse();
       Vector2s V1p = V1 + M12 * Vaux;
       Vector2s V2p = V2 - M21 * Vaux;
       */

    Matrix2s M1 = M.segment<2>(I1).asDiagonal();
    Matrix2s M2 = M.segment<2>(I2).asDiagonal();
    bool F1 = scene.isFixed(idx1);
    bool F2 = scene.isFixed(idx2);

    calculateVertexVertexScalarMInfResponse(I1, I2, V1, V2, Vaux, N, M1, M2, F1, F2, V, V1p, V2p);

#ifndef NDEBUG
    cout << "**** " << __FUNCTION__ \
        << "\n[" << g_Time << "]  PP collision: n:(" << n.transpose() << ") COR=" << getCOR() << " Vaux:(" << Vaux //Vaux.transpose()
        << ") V1:(" << V1.transpose() << "), V1p:(" << V1p.transpose() \
        << ") V2:(" << V2.transpose() << "), V2p:(" << V2p.transpose() << ")" \
        << ", M2:" << M2 << ", M1:" << M1 \
        << endl;
#endif
}

// Responds to a collision detected between a particle and an edge by applying
// an impulse to the velocities of each one.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge.
//   n:     The shortest vector between the particle and the edge.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleEdge(TwoDScene &scene, int vidx, int eidx, const Vector2s &n)
{
    // Your code goes here!
    pair<int, int> const &E = scene.getEdges()[eidx];

    int I1 = vidx << 1;
    int I2 = E.first << 1;
    int I3 = E.second << 1;

    VectorXs const &X = scene.getX();
    VectorXs const &M = scene.getM();
    VectorXs &V = scene.getV();

    Vector2s X1 = X.segment<2>(I1);
    Vector2s X2 = X.segment<2>(I2);
    Vector2s X3 = X.segment<2>(I3);

    Vector2s V1 = V.segment<2>(I1);
    Vector2s V2 = V.segment<2>(I2);
    Vector2s V3 = V.segment<2>(I3);

    Matrix2s M1 = M.segment<2>(I1).asDiagonal();
    Matrix2s M2 = M.segment<2>(I2).asDiagonal();
    Matrix2s M3 = M.segment<2>(I3).asDiagonal();
    Vector2s V1p, V2p, V3p, N;
    bool F1 = scene.isFixed(vidx);
    bool F2 = scene.isFixed(E.first);
    bool F3 = scene.isFixed(E.second);
    V1p = V2p = V3p = Vector2s::Zero();
    // alpha parameter (edge drop point absolute distance)
    Vector2s Ve;
    scalar sAlpha;

    findAlpha(X1, X2, X3, V1, V2, V3, N, Ve, sAlpha);

    N = n; // calculated from detection not now
    N.normalize();
    scalar C = (1.0 + getCOR());
    scalar Vaux = C * (Ve-V1).dot(N);

    calculateVertexEdgeScalarMInfResponse(I1, I2, I3, V1, V2, V3, Vaux, N, M1, M2, M3, F1, F2, F3, sAlpha, V, V1p, V2p, V3p);

#ifndef NDEBUG
    cout << "**** " << __FUNCTION__ \
        << "\n[" << g_Time << "]  PP collision: n:(" << n.transpose() << ") COR=" << getCOR() << " Alpha:" << sAlpha \
        << " Vaux:(" << Vaux //Vaux.transpose()
        << ") V1:(" << V1.transpose() << "), V1p:(" << V1p.transpose() \
        << ") V2:(" << V2.transpose() << "), V2p:(" << V2p.transpose() << ")" \
        << ") V3:(" << V2.transpose() << "), V3p:(" << V3p.transpose() << ")" \
        << ", M1:" << M1 << ", M2:" << M2 << ", M3:" << M3 \
        << endl;
#endif
}


// Responds to a collision detected between a particle and a half-plane by 
// applying an impulse to the velocity of the particle.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the half-plane.
//   n:     The shortest vector between the particle and the half-plane.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleHalfplane(TwoDScene &scene, int vidx, int pidx, const Vector2s &n)
{
    // Your code goes here!
    VectorXs &V = scene.getV();
    Vector2s Xp = scene.getHalfplane(pidx).first;
    bool F1 = scene.isFixed(vidx);

    if (F1)
        return;

    int I1 = vidx << 1;
    scalar C = (1.0 + getCOR());
    Vector2s V1 = V.segment<2>(I1);
    Vector2s N = n;
    N.normalize(); // just in case

    Vector2s Vaux = C * (V1).dot(N) * N;
    Vector2s V1p = V1 - Vaux;

    V.segment<2>(I1) = V1p;

#ifndef NDEBUG
    cout << "**** " << __FUNCTION__ \
        << "\n[" << g_Time << "]  PP collision: n:(" << n.transpose() << ") COR=" << getCOR() << " Vaux:(" << Vaux \
        << ") V1:(" << V1.transpose() << "), V1p:(" << V1p.transpose() \
        << endl;
#endif

}
