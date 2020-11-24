#include "RigidBodyWindForce.h"

#include "../Ops.h"


static void computeEdgeWind(
    RigidBody &_RB,
    Vector2s const &_Xei,
    Vector2s const &_Edge,
    int _PNum,
    scalar const &_Beta,
    Vector2s const &_Wind,
    Vector2s &_EdgeN,
    scalar &_T,
    Vector2s &_F)
{
  scalar EdgeLen = _Edge.norm();
  Vector2s EdgeHat = _Edge; EdgeHat.normalize();
  Matrix2s R90; R90 << 0,-1,1,0;
  Vector2s EdgeN = R90 * EdgeHat; EdgeN.normalize();
  scalar P = EdgeLen / _PNum;
  Vector2s Xi = _Xei + P/2 * EdgeHat;

  scalar A = _Beta * P;

  for (int i=0; i < _PNum; ++i) {
    Vector2s Vi = _RB.computeWorldSpaceVelocity(Xi);
    Vector2s F = A * (_Wind - Vi).dot(EdgeN) * EdgeN;
    Vector2s Ri = (Xi - _RB.getX());
    Vector2s RiN = Ri; RiN.normalize();
    scalar T = cross2s(Ri, F); // ?? R90*RiN.dot(F) * Ri.norm() ??

    Xi += P * EdgeHat;
    _F += F;
    _T += T;
  };

  _EdgeN = EdgeN;
}



static void calculateRB(
    RigidBody &_RB,
    int _PNum,
    scalar const &_Beta,
    Vector2s const &_Wind)
{
  int Edges = _RB.getNumEdges();
  Vector2s F = Vector2s::Zero();
  Vector2s EdgeN = Vector2s::Zero();
  scalar T = 0.0;

  for (int i=0; i < Edges; ++i) {
    computeEdgeWind(
        _RB,
        _RB.getWorldSpaceVertex(i),
        _RB.computeWorldSpaceEdge(i),
        _PNum,
        _Beta,
        _Wind,
        EdgeN,
        T,
        F);
  }

  _RB.getForce() += F;
  _RB.getTorque() += T;
}




void RigidBodyWindForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  // Your code goes here!
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  //
  int RBNum = rbs.size();
  for (int i=0; i < RBNum; ++i) {
    calculateRB(rbs[i], m_num_quadrature_points, m_beta, m_wind);
  }
}
