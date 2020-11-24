#include "RigidBodyWindForce.h"



static void computeEdgeWind(
    Vector2s const &_Edge,
    scalar const &_Beta,
    Vector2s const &_Wind)
{
}



static void calculateRB(
    RigidBody &_RB,
    scalar const &_Beta,
    Vector2s const &_Wind)
{
  int Edges = _RB.getNumEdges();

  for (int i=0; i < Edges; ++i) {
    computeEdgeWind(_RB.computeWorldSpaceEdge(i), _Beta, _Wind);
  }
}




void RigidBodyWindForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  // Your code goes here!
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  //


  calculateRB(rbs[0], m_beta, m_wind);

}
