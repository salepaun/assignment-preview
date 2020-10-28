#include "SimpleCollisionHandler.h"
#include <iostream>
#include <set>
#include <vector>
// #include <lapacke.h>
#include <Eigen/Dense>

#include "MathDefs.h"



using namespace std;



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

    if(bTouch)
      bApproach = (V1 - V2).dot(n) > 0;
    else
      bApproach = false;

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << "\n  X1:" << X1.transpose() << ", X2:" << X2.transpose() \
      << ", D:" << D << ", R1:" << R1 << ", R2:" << R2 \
      << ": T:" << bTouch << ", A:" << bApproach << endl;
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
    scalar R = max(max(R2,R3),Re); //max3(R2,R3,Re);

    Vector2s A = X1 - X2;
    Vector2s B = X3 - X2;
    Vector2s EAxN = B;
    EAxN.normalize();

    Vector2s X1EAx;
    Vector2s Vx;

    // alpha parameter (edge drop point absolute distance)
    scalar BLen = B.norm();
    scalar sEAlpha = A.dot(B) / BLen;
    scalar sAlpha = sEAlpha / BLen;

    // 
    // Alpha constraints
    //
    if (sAlpha <= 0) {
      // away from or on X2 edge vertex
      X1EAx = (X1 - X2);
      Vx = V2;
    } else if (sAlpha >= 1) {
      // away from or on X3 edge vertex
      X1EAx = (X1 - X3);
      Vx = V3;
    } else {
      // calculating Normal and distance with 'drop' point
      X1EAx = sEAlpha * EAxN - A;
      Vx = V2 + sAlpha * (V3 - V2);
    };

    // output according to description (not normalized!)
    n = -X1EAx;

    scalar D = X1EAx.norm();
    X1EAx.normalize(); // need normalized for speed vectors

    bTouch = D < R1 + R;

    if(bTouch)
      bApproach = (V1 - Vx).dot(n) > 0;
    else
      bApproach = false;
    
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

    if(bTouch)
      bApproach = V1.dot(n) > 0;
    else
      bApproach = false;

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << ", X1:" << X1.transpose() << ", Xp:" << Xp.transpose() \
      << ", D:" << D << ", R1:" << R1 \
      << ": T:" << bTouch << ", A:" << bApproach << endl;
#endif

    return bTouch && bApproach;
    
    return false;
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
    scalar C = (1.0 + getCOR());
    Vector2s V1 = V.segment<2>(I1);
    Vector2s V2 = V.segment<2>(I2);
    Vector2s N = n;
    N.normalize();

    // n.normalize(); // just in case
    Vector2s Vaux = C * (V2-V1).dot(N) * N;

    Matrix2s M1 = M.segment<2>(I1).asDiagonal();
    Matrix2s M2 = M.segment<2>(I2).asDiagonal();
    Matrix2s M12 = M2 + M1;
    Matrix2s M12i = M12.inverse();
    Vector2s V1p = V1 + M2 * M12i * Vaux;
    Vector2s V2p = V2 - M1 * M12i * Vaux;

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

    V.segment<2>(I1) = V1p;
    V.segment<2>(I2) = V2p;

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << "\n  PP collision: n:(" << n.transpose() << ") COR=" << getCOR() << " Vaux:(" << Vaux \
      << ") V1:(" << V1.transpose() << "), V1p:(" << V1p.transpose() \
      << ") V2:(" << V2.transpose() << "), V2p:(" << V2p.transpose() << ")" \
      << ", M2:" << M2 << ", M1:" << M1 << ", M12i:" << M12i \
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
    const VectorXs &M = scene.getM();
    
    int eidx1 = scene.getEdges()[eidx].first;
    int eidx2 = scene.getEdges()[eidx].second;
    
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*eidx1);
    VectorXs x3 = scene.getX().segment<2>(2*eidx2);
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*eidx1);
    VectorXs v3 = scene.getV().segment<2>(2*eidx2);
    
    // Your code goes here!
    
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
    VectorXs nhat = n;
    
    // Your code goes here!
    
}
