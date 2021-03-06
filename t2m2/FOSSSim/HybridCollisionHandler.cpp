#include "HybridCollisionHandler.h"
#include "ContinuousTimeUtilities.h"
#include <iostream>
#include <set>
#include <algorithm>
#include "HybridCollisionComparison.h"

#include "MathDefs.h"
#include "Ops.h"



using namespace std;



/**********************************************************************************
 * Local not exported methods.
 */


inline static double cross(Vector2s _A, Vector2s _B)
{
    return _A[0] * _B[1] - _A[1] * _B[0]; // or through M.det
}




static void calculateImpactZone(
        TwoDScene const &_Scene,
        ImpactZone const &_Zone,
        VectorXs const &_Qs,
        VectorXs &_dQ,
        Vector2s &_ZoneCM,
        double &_ZoneRot,
        Vector2s &_ZoneDeltaCM)
{
    VectorXs const &M = _Scene.getM();

    Vector2s XMsum = Vector2s::Zero();
    Vector2s dQMsum = Vector2s::Zero();

    double Lsum = 0;
    double Isum = 0;
    double Msum = 0;

    set<int>::const_iterator i = _Zone.m_verts.begin();

    // mass calculations on scalar m1 (following answer to a previous assignment)
    for (;i != _Zone.m_verts.end(); ++i) { // for (auto i: _Zone.m_verts) {
        int j = (*i) <<1;
        double m = M.segment<2>(j)[0];
        Msum += m;
        XMsum += m * _Qs.segment<2>(j);
        dQMsum += m * _dQ.segment<2>(j);
    };

    _ZoneCM = XMsum / Msum;
    _ZoneDeltaCM = dQMsum / Msum;

    for (i=_Zone.m_verts.begin(); i != _Zone.m_verts.end(); ++i) { // for (auto i: _Zone.m_verts) {
        int j = (*i)<<1;
        double m = M.segment<2>(j)[0];
        Vector2s xi = _Qs.segment<2>(j);
        Vector2s vi = _dQ.segment<2>(j);
        Lsum += m * cross((xi - _ZoneCM), (vi - _ZoneDeltaCM));
        Isum += m * (xi - _ZoneCM).squaredNorm();
    };

    _ZoneRot = Lsum / Isum;
}



static void applyImpactZone(
        double _Dt,
        ImpactZone const &_Zone,
        VectorXs const &_Qs,
        Vector2s const &_ZoneCM,
        double const &_ZoneRot,
        Vector2s const &_ZoneDeltaCM,
        VectorXs &_Qe,
        VectorXs &_Ve)
{
    set<int>::const_iterator i = _Zone.m_verts.begin();
    for (;i != _Zone.m_verts.end(); ++i) { // for (auto i: _Zone.m_verts) {
        int j = (*i)<<1;
        Vector2s Xi = _Qs.segment<2>(j);

        Vector2s Xir = (Xi - _ZoneCM);
        Vector2s Xiro(-Xir[1], Xir[0]);
        Vector2s Xie = _ZoneCM + _ZoneDeltaCM + cos(_ZoneRot)*(Xir) + sin(_ZoneRot)*(Xiro);

        _Qe.segment<2>(j) = Xie;
        _Ve.segment<2>(j) = (Xie - Xi) / _Dt;
    };
}







////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////        Impact Zone Utilities
////
////        You can use them (but please do make sure you understand what they do), or 
////        implement your own versions of whatever you need
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

bool intersects(const ImpactZone &z1, const ImpactZone &z2)
{
    std::set<int> temp;
    set_intersection(z1.m_verts.begin(), z1.m_verts.end(), z2.m_verts.begin(), z2.m_verts.end(), inserter(temp, temp.end()));
    return temp.size() > 0;
}

ImpactZone mergeZones(const ImpactZone &z1, const ImpactZone &z2)
{
    std::set<int> combinedverts;
    set_union(z1.m_verts.begin(), z1.m_verts.end(), z2.m_verts.begin(), z2.m_verts.end(), inserter(combinedverts, combinedverts.end()));
    return ImpactZone(combinedverts, z1.m_halfplane || z2.m_halfplane);
}


void mergeAllZones(ImpactZones &zones)
{
    ImpactZones result;
    
    ImpactZones *src = &zones;
    ImpactZones *dst = &result;
    do
    {
        dst->clear();
        for(int i=0; i<(int)src->size(); i++)
        {
            bool merged = false;
            for(int j=0; j<(int)dst->size(); j++)
            {
                if(intersects((*dst)[j], (*src)[i]))
                {
                    ImpactZone newzone = mergeZones((*dst)[j], (*src)[i]);
                    (*dst)[j] = newzone;
                    merged = true;
                    
                    break;
                }
            }
            if(!merged)
            {
                dst->push_back((*src)[i]);
            }
        }
        std::swap(src, dst);
    }
    while(src->size() < dst->size());
    
    zones = *dst;
}

void growImpactZones(const TwoDScene &scene, ImpactZones &zones, const std::vector<CollisionInfo> &impulses)
{
    for(int i=0; i<(int)impulses.size(); i++)
    {
        switch(impulses[i].m_type)
        {
            case CollisionInfo::PP:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                verts.insert(impulses[i].m_idx2);
                zones.push_back(ImpactZone(verts, false));
                break;
            }
            case CollisionInfo::PE:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                verts.insert(scene.getEdge(impulses[i].m_idx2).first);
                verts.insert(scene.getEdge(impulses[i].m_idx2).second);
                zones.push_back(ImpactZone(verts, false));
                break;
            }
            case CollisionInfo::PH:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                zones.push_back(ImpactZone(verts, true));
                break;
            }
        }
    }
    mergeAllZones(zones);
}

bool zonesEqual(const ImpactZones &zones1, const ImpactZones &zones2)
{
    if(zones1.size() != zones2.size())
        return false;
    
    for(int i=0; i<(int)zones1.size(); i++)
    {
        bool found = false;
        for(int j=0; j<(int)zones2.size(); j++)
        {
            if(zones1[i] == zones2[j])
            {
                found = true;
                
                break;
            }
        }
        if(!found)
            return false;
    }
    return true;
}





////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////        Student Code
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////



// Iteratively performs collision detection and interative impulse response until either there are no more detected collisions, or the maximum number of
// iterations has been reached. See the assignment instructions for more details.
// The maximum number of iterations is stored in the member variable m_maxiters.
// Inputs:
//   scene:   The simulation scene. Get masses, radii, edge endpoint indices, etc. from here. Do *NOT* get any positions or velocities from here.
//   qs:      The positions of the particles at the start of the time step.
//   qe:      The predicted end-of-time-step positions.
//   qdote:   The predicted end-of-time-step velocities.
//   dt:      The time step size.
// Outputs:
//   qefinal:    The collision-free end-of-time-step positions (if no new collisions are detected), or the last set of predicted end-of-time-step positions
//               (if maximum number of iterations reached).
//   qdotefinal: Same as qefinal, but for velocities.
//   Returns true if the algorithm found a collision-free state. Returns false if the maximum number of iterations was reached without finding a collision-
//   free state.
// Possibly useful functions: detectCollisions, applyImpulses.
bool HybridCollisionHandler::applyIterativeImpulses(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, const VectorXs &qdote, double dt, VectorXs &qefinal, VectorXs &qdotefinal)
{
    // Your code goes here!
    bool bRet = false;
    int i = 0;

    vector<CollisionInfo> Collisions;

    VectorXs qe_tmp;
    VectorXs qdote_tmp;

    qefinal = qe;
    qdotefinal = qdote;

    for (i=0; i < m_maxiters; ++i) {
        Collisions = detectCollisions(scene, qs, qefinal);
        if (Collisions.size()) {
            applyImpulses(scene, Collisions, qs, qefinal, qdote, dt, qe_tmp, qdote_tmp);
            qefinal = qe_tmp;
            qdotefinal = qdote_tmp;
        } else {
            bRet = true;
            break;
        };
    };

#ifndef NDEBUG
    cout << __FUNCTION__
        << ", Iterations:=" << i << ", max=" << m_maxiters << " ret:" << bRet
        << " Collisions num=" << Collisions.size()
        << endl;
#endif

    return bRet;
}


// Resolves any remaining collisions in a simulation time step by setting the velocities of all particles involved in a way that guarantees
// that the distance between particles in an impact zone does not change.
// Inputs:
//   scene:   The simulation scene, from which the masses of the particles, current (colliding) positions, and whether or not a given particle is fixed,
//             can be retrieved.
//   qs:      The positions of the particles at the start of the time step.
//   qe:      The predicted end-of-timestep positions of the particles.
//   qdote:   The precicted end-of-timestep velocities of the particles.
//   zone:    Information about the impact zone of colliding particles. zone.m_verts is an std::set of particle indices; each particle in this set
//            is part of the impact zone and needs to have its position and velocity changed. Whether or not a half-plane is part of the impact zone
//            can be checked by looking at zone.m_halfplane.
//   dt:      The time step.
// Outputs:
//   qe:      The end-of-timestep position of the particles, assuming rigid motion as in writeup section 4.5
//   qdote:   The end-of-timestep velocity of the particles, assuming rigid motion as in writeup section 4.5
void HybridCollisionHandler::performFailsafe(const TwoDScene &scene, const VectorXs &qs, const ImpactZone &zone, double dt, VectorXs &qe, VectorXs &qdote)
{
    //
    // What you need to implement here: (same as writeup section 4.5)
    //
    // 1. Treat the particles as if they were part of a rigid body; that is, treat them as if 
    //      we connected them with rigid beams at the start of the time step.
    // 2. Step the rigid body forward in time to the end of the time step.
    // 3. Set each particle’s modified end-of-time-step position qm to the position dictated 
    //      by the motion of the rigid body.
    // 4. Also set the particle’s modified end-of-time-step velocity to (qm − qs) / h, where 
    //      h is the length of the time step.
    //
    
    // Don't forget to handle fixed objects properly as in writeup section 4.5.1    
    
    // Your code goes here!

    //SetVertedFixed SetVertexFixedOp(qs, qe, qdote);
    // IsVertexFixed IsVertexFixedPred(scene);
    // bool bFixedZone = zone.m_halfplane || std::any_of(zone.m_verts.begin(), zone.m_verts.end(), IsVertexFixedPred);
    bool bFixedZone = zone.m_halfplane;

    VectorXs dQ = qe - qs;
    Vector2s ZoneCM = Vector2s::Zero();
    double ZoneRot = 0;
    Vector2s ZoneDeltaCM = Vector2s::Zero();

    for (set<int>::const_iterator i=zone.m_verts.begin(); i!=zone.m_verts.end(); ++i) {
        if (scene.isFixed(*i)) {
            bFixedZone = true;
            break;
        };
    };

    // Impact zone involved with fixed objects as per this milestone - just freeze
    if (bFixedZone) {
        for (set<int>::const_iterator i=zone.m_verts.begin(); i != zone.m_verts.end(); ++i) { // : zone.m_verts) {
            int j = (*i)<<1;
            qe.segment<2>(j) = qs.segment<2>(j);
            qdote.segment<2>(j).setZero();
        };
    } else {
        calculateImpactZone(scene, zone, qs, dQ, ZoneCM, ZoneRot, ZoneDeltaCM);
        applyImpactZone(dt, zone, qs, ZoneCM, ZoneRot, ZoneDeltaCM, qe, qdote);
    };

#ifndef NDEBUG
    std::cout << __FUNCTION__
        << ", Zone Vnum=" << zone.m_verts.size() << ", HalfPlane:" << zone.m_halfplane
        << ", Fixed:" << bFixedZone
        << ", ZoneCM:" << ZoneCM.transpose()
        << ", ZoneRot:" << ZoneRot
        << ", ZoneDeltaCM:" << ZoneDeltaCM.transpose()
        << ", dt=" << dt
        << std::endl;
#endif

}


// Performs iterative geometric collision response until collision-free end-of-time-step positions and velocities are found. See the assignment
// instructions for details.
// Inputs:
//   scene:   The simulation scene. Get masses, radii, etc. from here. Do *NOT* get any positions or velocities from here.
//   qs:      The start-of-time-step positions.
//   qe:      The predicted end-of-time-step positions.
//   qdote:   The predicted end-of-time-step velocities.
//   dt:      The time step size.
// Outputs:
//   qm:    The final, collision-free end-of-time-step positions. (qm in the assignment instructions.)
//   qdotm: Same as qm, but for velocities.
// Possibly useful functions: detectCollisions, performFailsafe. You may find it helpful to write other helper functions for manipulating (merging,
// growing, etc) impact zones.
void HybridCollisionHandler::applyGeometricCollisionHandling(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, const VectorXs &qdote, double dt, VectorXs &qm, VectorXs &qdotm)
{
    ImpactZones Z;
    ImpactZones Zprime;
    
    //
    // What you need to implement here: (same as writeup section 4.6)
    //
    // 1. Perform continuous-time collision detection using positions qs and qe.
    // 2. Initialize qm = qe and qdotm = qdote.
    // 3. Construct a list of disjoint impact zones Z from the detected collisions.
    // 4. For each impact zone in Z, apply geometric collision response (by calling 
    //      HybrdiCollisionHandler::performFailsafe, using positions qs and qm, and
    //      modifying qm and qdotm for the vertices in those zones.
    // 5. Perform continuous-time collision detection using positions qs and qm. 
    // 6. Construct a new list of impact zones Z′ consisting of all impact zones in Z, 
    //      plus one zone for each detected collision.
    // 7. Merge the zones in Z' to get disjoint impact zones.
    // 8. If Z and Z' are equal, the algorithm is done, and qm and qdotm are the new, 
    //      collision-free end-of-time-step positions. Z and Z' are equal if they 
    //      contain exactly the same impact zones; impact zones are the same if they 
    //      contain the same particles and they both involve, or both don’t involve, 
    //      a half-plane. If Z != Z', go to step 9.
    // 9. Set Z=Z' and goto step4.
    //
    
    // Your code goes here!
    vector<CollisionInfo> Collisions;

    int i = 0;
    bool bContinue = true;
    Collisions = detectCollisions(scene, qs, qe);
    qm = qe;
    qdotm = qdote;
    growImpactZones(scene, Z, Collisions);

    for (i=0; bContinue; ++i) {
        for (vector<ImpactZone>::iterator i=Z.begin(); i!=Z.end(); ++i) { // auto Zone: Z) {
            ImpactZone &Zone = (*i);
            performFailsafe(scene, qs, Zone, dt, qm, qdotm);
        };
        Collisions = detectCollisions(scene, qs, qm);
        Zprime = Z;
        growImpactZones(scene, Zprime, Collisions);

        if (zonesEqual(Z, Zprime)) {
            bContinue = false;
            continue;
        }

        Z = Zprime;
    };

    /*
    do {
        Z = Zprime;
        Collisions = detectCollisions(scene, qs, qm);
        growImpactZones(scene, Zprime, Collisions);
        if (Collisions.size()) {
            for (auto Zone: Zprime) {
                performFailsafe(scene, qs, Zone, dt, qm, qdotm);
            };
        };
    } while (!zonesEqual(Z, Zprime));
    */

#ifndef NDEBUG
    cout << __FUNCTION__
        << " Collisions num=" << Collisions.size()
        << " Zones num=" << Z.size()
        << " dt=" << dt
        << endl;
#endif

}

