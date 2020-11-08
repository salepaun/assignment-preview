#include "ContinuousTimeCollisionHandler.h"
#include <iostream>
#include "ContinuousTimeUtilities.h"

// BEGIN STUDENT CODE //

#include <iomanip>

#include "Ops.h"


using namespace std;


// Given the start position (oldpos) and end position (scene.getX) of two
// particles, and assuming the particles moved in a straight line between the
// two positions, determines whether the two particles were overlapping and
// approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the two particles
// at the time of collision, and sets t to the time (0 = start position, 
// 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. The new positions and radii of the 
//           particles can be obtained from here.
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   idx1:   The index of the first particle. (ie, the degrees of freedom
//           corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//           scene.getX()).
//   idx2:   The index of the second particle.
// Outputs:
//   n:    The vector between the two particles at the time of collision.
//   time: The time (scaled to [0,1]) when the two particles collide.
//   Returns true if the two particles overlap and are approaching at some point
//   during the motion.
bool ContinuousTimeCollisionHandler::detectParticleParticle(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int idx1, int idx2, Vector2s &n, double &time)
{
    VectorXs dx = qe-qs;
    
    Vector2s X1 = qs.segment<2>(idx1<<1);
    Vector2s X2 = qs.segment<2>(idx2<<1);

    Vector2s dX1 = dx.segment<2>(idx1<<1);
    Vector2s dX2 = dx.segment<2>(idx2<<1);
    
    double R1 = scene.getRadius(idx1);
    double R2 = scene.getRadius(idx2);

    vector<Polynomial> polynomials;
    vector<double> position_polynomial;
    vector<double> velocity_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients
  
    position_polynomial.push_back(-(dX2-dX1).dot(dX2-dX1));
    position_polynomial.push_back(-2*(X2-X1).dot(dX2-dX1));
    position_polynomial.push_back((R1+R2)*(R1+R2) - (X2-X1).dot(X2-X1));
    
    velocity_polynomial.push_back((dX1-dX2).dot(dX2-dX1));
    velocity_polynomial.push_back((dX1-dX2).dot(X2-X1));
    /*
    position_polynomial.push_back(-1);
    position_polynomial.push_back(1);
    position_polynomial.push_back(0);
    velocity_polynomial.push_back(1);
    velocity_polynomial.push_back(-0.5);
    */

    // Do not change the order of the polynomials here, or your program will fail the oracle
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(velocity_polynomial));

    // Your implementation here should compute n, and examine time to decide the return value
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    n = (X2 + time * dX2) - (X1 + time * dX1);
    bool bCollision = time >= 0 && time <= 1;

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(5)
      << (bCollision ? "COLLISION" : " no collision")
      << ", X1=" << X1.transpose() << ", dX1:" << dX1.transpose()
      << ", X2=" << X2.transpose() << ", dX2:" << dX2.transpose()
      << ", N=" << n.transpose() << ", T=" << time << endl;
#endif

      return bCollision;
}


// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and an edge, and assuming the particle and edge endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the particle and the
// edge at the time of collision, and sets t to the time (0 = start position, 
// 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. 
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   vidx:   The index of the particle.
//   eidx:   The index of the edge.
// Outputs:
//   n:    The shortest vector between the particle and edge at the time of 
//         collision.
//   time: The time (scaled to [0,1]) when the two objects collide.
//   Returns true if the particle and edge overlap and are approaching at
//   some point during the motion.
// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and an edge, and assuming the particle and edge endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
bool ContinuousTimeCollisionHandler::detectParticleEdge(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int vidx, int eidx, Vector2s &n, double &time)
{
    VectorXs dx = qe - qs;
    
    VectorXs x1 = qs.segment<2>(vidx<<1);
    VectorXs x2 = qs.segment<2>(scene.getEdge(eidx).first<<1);
    VectorXs x3 = qs.segment<2>(scene.getEdge(eidx).second<<1);
    
    VectorXs dx1 = dx.segment<2>(vidx<<1);
    VectorXs dx2 = dx.segment<2>(scene.getEdge(eidx).first<<1);
    VectorXs dx3 = dx.segment<2>(scene.getEdge(eidx).second<<1);

    double r1 = scene.getRadius(vidx);
    double r2 = scene.getEdgeRadii()[eidx];

    std::vector<double> position_polynomial;
    std::vector<double> alpha_greater_than_zero_polynomial;
    std::vector<double> alpha_less_than_one_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients

    Vector2s X1 = x1.segment<2>(0);
    Vector2s X2 = x2.segment<2>(0);
    Vector2s X3 = x3.segment<2>(0);

    Vector2s dX1 = dx1.segment<2>(0);
    Vector2s dX2 = dx2.segment<2>(0);
    Vector2s dX3 = dx3.segment<2>(0);

    double rr2 = (r1+r2)*(r1+r2);

    Vector2s dX12 = dX1-dX2;
    Vector2s dX32 = dX3-dX2;
    Vector2s X12 = X1-X2;
    Vector2s X32 = X3-X2;

    Vector2s dX13 = dX1-dX3;
    Vector2s dX23 = dX2-dX3;
    Vector2s X13 = X1-X3;
    Vector2s X23 = X2-X3;
    
    double dX12_dX32 = dX12.dot(dX32);
    double dX12_dX12 = dX12.dot(dX12);
    double dX32_dX32 = dX32.dot(dX32);
    double X12_dX12 = X12.dot(dX12);
    double X32_dX32 = X32.dot(dX32);
    double X12_dX32 = X12.dot(dX32);
    double dX12_X32 = dX12.dot(X32);
    double X32_X32 = X32.dot(X32);
    double X12_X12 = X12.dot(X12);
    double X12_X32 = X12.dot(X32);
    double X32_dX12 = X32.dot(dX12); //dX12_X32;
    double A1 = (X12_dX32 + dX12_X32);
 
    // t4
    position_polynomial.push_back(
        dX12_dX32*dX12_dX32 - dX12_dX12*dX32_dX32);
    // t3
    position_polynomial.push_back(-2*(
        dX32_dX32*X12_dX12 + dX12_dX12*X32_dX32 - dX12_dX32*X12_dX32 + dX12_X32));
    // t2
    position_polynomial.push_back(
        2 * rr2 * dX32_dX32 
        - X32_X32 * dX12_dX12
        - X12_X12 * dX32_dX32
        - 4 * X12_dX12 * X32_dX32
        + 2 * X12_X32 * dX12_dX32
        + A1*A1);
    // t1
    position_polynomial.push_back(
        2 * (rr2 * X32_dX32
        - X32_X32 * X12_dX12
        - X12_X12 * X32_dX32
        + X12_X32 * (X12_dX32 + dX12_X32)));
    // t0
    position_polynomial.push_back(
        rr2 * X32_X32
        - X32_X32 * X12_X12
        + X12_X32 * X12_X32);


    alpha_greater_than_zero_polynomial.push_back(
        (dX12.dot(dX32)));
    alpha_greater_than_zero_polynomial.push_back(
        (X12.dot(dX32) + dX12.dot(X32)));
    alpha_greater_than_zero_polynomial.push_back(
        X12.dot(X32));

    alpha_less_than_one_polynomial.push_back(
        dX13.dot(dX23));
    alpha_less_than_one_polynomial.push_back(
        X13.dot(dX23)+dX13.dot(X23));
    alpha_less_than_one_polynomial.push_back(
        X13.dot(X23));

    // Here's the quintic velocity polynomial:
    std::vector<double> velcity_polynomial;
    {
        double a = (x3-x2).dot(x3-x2);
        double b = (x3-x2).dot(dx3-dx2);
        double c = (dx3-dx2).dot(dx3-dx2);
        double d = (dx2-dx1).dot(dx2-dx1);
        double e = (dx2-dx1).dot(x2-x1);
        double f = (x1-x2).dot(x3-x2);
        double g = (x1-x2).dot(dx3-dx2) + (dx1-dx2).dot(x3-x2);
        double h = (dx1-dx2).dot(dx3-dx2);
        double i = (dx3-dx2).dot(x2-x1) + (dx2-dx1).dot(x3-x2);
        double j = (dx3-dx2).dot(dx2-dx1);
        double k = a*f;
        double l = a*g+2*b*f;
        double m = a*h+2*b*g+c*f;
        double n = c*g+2*b*h;
        double o = c*h;
        double p = (dx3-dx2).dot(x3-x2);
        double q = (dx3-dx2).dot(dx3-dx2);
        
        velcity_polynomial.push_back( -h*h*q - c*c*d - 2*o*j );
        velcity_polynomial.push_back( -h*h*p - 2*g*h*q - 4*b*c*d - c*c*e - o*i - 2*n*j );
        velcity_polynomial.push_back( -2*g*h*p - 2*f*g*q - g*g*q - 2*a*c*d - 4*b*b*d - 4*b*c*e - n*i - 2*m*j );
        velcity_polynomial.push_back( -2*f*h*p - g*g*p - 2*f*g*q - 4*a*b*d - 2*a*c*e - 4*b*b*e - m*i - 2*l*j );
        velcity_polynomial.push_back( -2*f*g*p - f*f*q - a*a*d - 4*a*b*e - l*i - 2*k*j );
        velcity_polynomial.push_back( -f*f*p - a*a*e - k*i );
    }

    // Do not change the order of the polynomials here, or your program will fail the oracle
    std::vector<Polynomial> polynomials;
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(alpha_greater_than_zero_polynomial));
    polynomials.push_back(Polynomial(alpha_less_than_one_polynomial));
    polynomials.push_back(Polynomial(velcity_polynomial));
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    
    // Your implementation here should compute n, and examine time to decide the return value
    Vector2s X1t = X1 + dX1*time;
    Vector2s X2t = X2 + dX2*time;
    Vector2s X3t = x3 + dX3*time;
    Vector2s V0 = Vector2s::Zero();
    Vector2s Ve = Vector2s::Zero();
    double sAlpha = 0; //(X1t-X2t).dot(X3t-X21) / (X3t-X2t).squaredNorm();
    findAlpha(X1t, X2t, X3t, V0, V0, V0, n, Ve, sAlpha);
    bool bCollision = time >= 0 && time <= 1;

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(5)
      << (bCollision ? "COLLISION" : " no collision")
      << ", X1=" << X1.transpose() << ", dX1:" << dX1.transpose()
      << ", X2=" << X2.transpose() << ", dX2:" << dX2.transpose()
      << ", X3=" << X3.transpose() << ", dX3:" << dX3.transpose()
      << ", N=" << n.transpose() << ", T=" << time << endl;
#endif

      return bCollision;
}

// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and a half-plane, and assuming the particle endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the particle and the
// half-plane at the time of collision, and sets t to the time (0 = start 
// position, 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. 
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   vidx:   The index of the particle.
//   eidx:   The index of the half-plane. The vectors (px, py) and (nx, ny) can
//           be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n:    The shortest vector between the particle and half-plane at the time 
//         of collision.
//   time: The time (scaled to [0,1]) when the two objects collide.
//   Returns true if the particle and half-plane overlap and are approaching at
//   some point during the motion.
// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and a half-plane, and assuming the particle endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
bool ContinuousTimeCollisionHandler::detectParticleHalfplane(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int vidx, int pidx, Vector2s &n, double &time)
{
    VectorXs dx = qe - qs;
    
    VectorXs x1 = qs.segment<2>(vidx<<1);
    VectorXs dx1 = dx.segment<2>(vidx<<1);
    
    VectorXs xp = scene.getHalfplane(pidx).first;
    VectorXs np = scene.getHalfplane(pidx).second;
    
    double r = scene.getRadius(vidx);
 
    std::vector<double> position_polynomial;
    std::vector<double> velocity_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients
    Vector2s X1 = x1.segment<2>(0);
    Vector2s dX1 = dx1.segment<2>(0);
    Vector2s Xp = xp.segment<2>(0);
    Vector2s Np = np.segment<2>(0);
    Vector2s Xp1 = Xp - X1;

    double Xp1Np = Xp1.dot(Np);
    double dX1Np = dX1.dot(Np);
    double r2 = r*r;

    position_polynomial.push_back(
        -dX1Np*dX1Np);
    position_polynomial.push_back(
        2*Xp1Np*dX1Np);
    position_polynomial.push_back(
        r2*Np.squaredNorm() - Xp1Np*Xp1Np);

    velocity_polynomial.push_back(
        -dX1Np*dX1Np);
    velocity_polynomial.push_back(
        (Xp-X1).dot(Np));

    // Do not change the order of the polynomials here, or your program will fail the oracle
    std::vector<Polynomial> polynomials;
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(velocity_polynomial));
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    
    // Your implementation here should compute n, and examine time to decide the return value
    Vector2s X1t = X1 + dX1*time;
    n = (Xp - X1t).dot(Np) * Np;

    bool bCollision = time >= 0 && time <= 1;

#ifndef NDEBUG
    cout << __FUNCTION__ << setprecision(5)
      << (bCollision ? "COLLISION" : " no collision")
      << ", Xp=" << Xp.transpose() << ", Np:" << Np.transpose()
      << ", X1=" << X1.transpose() << ", dX1:" << dX1.transpose() << ", X1t:" << X1t.transpose()
      << ", N=" << n.transpose() << ", T=" << time << endl;
#endif

    return bCollision;
}


