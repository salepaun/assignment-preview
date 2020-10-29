#include "PenaltyForce.h"
#include "TwoDScene.h"

#include "Ops.h"


using namespace std;



void PenaltyForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    // Feel free to implement if you feel like doing so.
}

// Adds the gradient of the penalty potential (-1 * force) for a pair of 
// particles to the total.
// Read the positions of the particles from the input variable x. Radii can
// be obtained from the member variable m_scene, the penalty force stiffness 
// from member variable m_k, and penalty force thickness from member variable
// m_thickness.
// Inputs:
//   x:    The positions of the particles in the scene. 
//   idx1: The index of the first particle, i.e. the position of this particle
//         is ( x[2*idx1], x[2*idx1+1] ).
//   idx2: The index of the second particle.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-particle
//          gradient to this total gradient.
void PenaltyForce::addParticleParticleGradEToTotal(const VectorXs &x, int idx1, int idx2, VectorXs &gradE)
{
    // Your code goes here!
    int I1 = idx1 << 1;
    int I2 = idx2 << 1;

    VectorXs const &X = x;
    
    vector<scalar> const & Radii = m_scene.getRadii();

    Vector2s X1 = X.segment<2>(I1);
    Vector2s X2 = X.segment<2>(I2);

    scalar R1 = Radii[idx1];
    scalar R2 = Radii[idx2];

    Vector2s Xd = X2 - X1;
    Vector2s N = Xd;

    N.normalize();

    scalar LenVx = Xd.norm();
    scalar D = R1 + R2 + m_thickness;

    if (LenVx <= D)
    {
      Matrix2s I = Matrix2s::Identity();
      MatrixXs GradN(2, 4);
      VectorXs GradE(4);

      GradN.block<2,2>(0,0) = -I;
      GradN.block<2,2>(0,2) = I;
      
      GradE = m_k * (LenVx - D) * GradN.transpose() * N;

      gradE.segment<2>(I1) += GradE.segment<2>(0);
      gradE.segment<2>(I2) += GradE.segment<2>(2);

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << " VxLen=" << LenVx << ", D=" << D \
      << ", GradE:" << GradE.transpose() \
      << endl;
#endif
    };

}

// Adds the gradient of the penalty potential (-1 * force) for a particle-edge
// pair to the total.
// Read the positions of the particle and edge endpoints from the input
// variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   eidx: The index of the edge, i.e. the indices of the particle making up the
//         endpoints of the edge are given by m_scene.getEdge(eidx).first and 
//         m_scene.getEdges(eidx).second.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-edge
//          gradient to this total gradient.
void PenaltyForce::addParticleEdgeGradEToTotal(const VectorXs &x, int vidx, int eidx, VectorXs &gradE)
{
    // Your code goes here!
    pair<int, int> const &E = m_scene.getEdges()[eidx];

    int I1 = vidx << 1;
    int I2 = E.first << 1;
    int I3 = E.second << 1;

    VectorXs const &X = x;
    
    vector<scalar> const & Radii = m_scene.getRadii();

    Vector2s X1 = X.segment<2>(I1);
    Vector2s X2 = X.segment<2>(I2);
    Vector2s X3 = X.segment<2>(I3);

    scalar R1 = Radii[vidx];
    scalar R2 = Radii[E.first];
    scalar R3 = Radii[E.second];
    scalar Re = m_scene.getEdgeRadii()[eidx];
    scalar R = Re; // seems like oracle uses edge, max(max(R2,R3),Re); //max3(R2,R3,Re);


    // alpha parameter (edge drop point absolute distance)
    Vector2s V1,V2,V3,Vx, N;
    scalar sAlpha;

    V1 = V2 = V3 = Vector2s::Zero();
    findAlpha(X1, X2, X3, V1, V2, V3, N, Vx, sAlpha, false);

    scalar D = N.norm();
    scalar Df = R1 + R + m_thickness;

    N.normalize();

    if (D <= Df)
    {
      Matrix2s I = Matrix2s::Identity();
      MatrixXs GradN(2, 6);
      VectorXs GradE(4);

      GradN.block<2,2>(0,0) = -I;
      GradN.block<2,2>(0,2) = (1-sAlpha) * I;
      GradN.block<2,2>(0,4) = sAlpha * I;
      
      GradE = m_k * (D - Df) * GradN.transpose() * N;

      gradE.segment<2>(I1) += GradE.segment<2>(0);
      gradE.segment<2>(I2) += GradE.segment<2>(2);
      gradE.segment<2>(I3) += GradE.segment<2>(4);

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << " D=" << D << ", Df=" << Df \
      << ", GradE:" << GradE.transpose() \
      << endl;
#endif
    };
    
}

// Adds the gradient of the penalty potential (-1 * force) for a particle-
// half-plane pair to the total.
// Read the positions of the particle from the input variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   pidx: The index of the half-plane, i.e. the position and normal vectors
//         for the half-plane can be retrieved by calling
//         m_scene.getHalfplane(pidx).
// Outputs:
//   gradE: The total gradient of the penalty force. *ADD* the particle-
//          half-plane gradient to this total gradient.
void PenaltyForce::addParticleHalfplaneGradEToTotal(const VectorXs &x, int vidx, int pidx, VectorXs &gradE)
{
    // Your code goes here!
    int I1 = vidx << 1;

    VectorXs const &X = x;

    VectorXs Xp = m_scene.getHalfplane(pidx).first;
    VectorXs Np = m_scene.getHalfplane(pidx).second;
    
    Vector2s X1 = X.segment<2>(I1);
    scalar R1 = m_scene.getRadii()[vidx];

    Np.normalize(); // just in case
    Vector2s N = (Xp - X1).dot(Np) * Np;

    scalar D = N.norm();
    scalar Df = R1 + m_thickness;

    N.normalize();

    if (D <= Df)
    {
      Matrix2s GradN;
      Vector2s GradE;

      GradN = -Np * Np.transpose() / Np.squaredNorm();
      GradE = m_k * (D - Df) * GradN.transpose() * N;
      gradE.segment<2>(I1) += GradE;

#ifndef NDEBUG
    cout << __FUNCTION__ \
      << " D=" << D << ", Df=" << Df \
      << ", GradE:" << GradE.transpose() \
      << endl;
#endif
    };
}
