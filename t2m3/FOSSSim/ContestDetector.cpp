/**
 */

#include <iostream>
#include <set>
#include <utility>
#include <deque>
#include <algorithm>
#include <iterator>
#include <array>

#ifndef NDEBUG
#include <assert.h>
#endif

#include "TwoDScene.h"
#include "ContestDetector.h"
#include "Ops.h"


using namespace std;



#define MY_DEBUG


/**********************************************************************************
 * ################################################################################
 * Local module declarations.
 */


typedef pair<double, double> T_MinMax;
typedef pair<int, int> T_IntersIds;
typedef set<T_IntersIds> T_IntersSet;


enum AlgoEnum {
  AllInAlgo = 0,
  EulerAlgo,
  LagrAlgo
};


/**********************************************************************************
 * ################################################################################
 * Local module variables.
 */


int g_MaxCoutNum = 3;



/**********************************************************************************
 * ################################################################################
 * Box.
 */



class Box {

  public:

    Box() {};
    Box(int _Id) : aId(_Id) {};

    inline int getId() const { return aId; };

    inline double const & xmin() const { return aBorder[0]; };
    inline double const & xmax() const { return aBorder[1]; };
    inline double const & ymin() const { return aBorder[2]; };
    inline double const & ymax() const { return aBorder[3]; };

    inline bool intersectX(Box const &_B) {
      return ((xmin() > _B.xmin()) && (xmin() < _B.xmax()))
        || ((xmax() > _B.xmin()) && (xmax() < _B.xmax()));
    };
    inline bool intersectY(Box const &_B) {
      return (ymin() > _B.ymin() && ymin() < _B.ymax())
        || (ymax() > _B.ymin() && ymax() < _B.ymax());
    };
    inline bool intersect(Box const &_B) {
      return intersectX(_B) && intersectY(_B);
    };

    friend ostream & operator << (ostream &_, Box const &_o);

  public:

    int aId;
    array<double, 4> aBorder;
};



ostream & operator << (ostream &_s, Box const &_o) {
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
  _s << "{" << _o.aId
    << " (nX:" << _o.aBorder[0]
    << ", mX:" << _o.aBorder[1]
    << ", nY:" << _o.aBorder[2]
    << ", mY:" << _o.aBorder[3]
    << ")";
  return _s;
}


/**********************************************************************************
 * ################################################################################
 * Boxed Object.
 */

/**
*/
class BoxedObj : public Box {

  public:

    BoxedObj() : Box() {};
    BoxedObj(int _Id) : Box(_Id) {};
    BoxedObj(int _Id, Vector2s const &_X, double const &_R) : Box(_Id), aX(_X[0]), aY(_X[1]) {
      // !!! aAvgSize += _R;
      updateR(_R);
    };

    inline bool changedX(Vector2s const &_X) {
      return _X[0] != aX || _X[1] != aY;
    };

    inline void updateR(double const &_R) {
      aBorder[0] = aX-_R;
      aBorder[1] = aX+_R;
      aBorder[2] = aY-_R;
      aBorder[3] = aY+_R;
    };

    inline void update(Vector2s const &_X, double const &_R) {
      aX = _X[0]; aY = _X[1];
      updateR(_R);
    };


    friend ostream & operator << (ostream &_, BoxedObj const &_o);


  private:

    double aX;
    double aY;
};



ostream & operator << (ostream &_s, BoxedObj const &_o) {
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
  _s << Box(_o)
    << " (" << _o.aX << "," << _o.aY
    << ")";
  return _s;
}





/**********************************************************************************
 * ################################################################################
 * Boxed Edge.
 */

/**
*/
class BoxedEdge : public Box {

  public:

    typedef pair<BoxedObj *, BoxedObj *> T_VrtxPair;

    BoxedEdge() : Box(), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id) : Box(_Id), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id, BoxedObj &_A, BoxedObj &_B, double const &_R) : 
      Box(_Id), aVrtx(&_A, &_B) {
        updateR(_R);
      };


    inline bool initialized() const {
      return aVrtx.first && aVrtx.second;
    };

    inline void updateR(double const &_R) {
      if (initialized()) {
        aBorder[0] = (a()->xmin() < b()->xmin() ? a()->xmin() : b()->xmin()) - _R;
        aBorder[1] = (a()->xmax() < b()->xmax() ? a()->xmax() : b()->xmax()) + _R;
        aBorder[2] = (a()->ymin() < b()->xmin() ? a()->xmin() : b()->xmin()) - _R;
        aBorder[3] = (a()->ymax() < b()->xmax() ? a()->xmax() : b()->xmax()) + _R;
      };
    };

    inline void update(double const &_R) {
      updateR(_R);
    };


    friend ostream & operator << (ostream &_, BoxedEdge const &_o);


  private:

    inline BoxedObj const * a() const {
      return aVrtx.first;
    };
    inline BoxedObj const * b() const {
      return aVrtx.second;
    };


  private:

    T_VrtxPair aVrtx;
};



ostream & operator << (ostream &_s, BoxedEdge const &_o) {
  _s << Box(_o)
    << " (A:" << (_o.a() ? _o.a()->getId() : -1)
    << ",B:" << (_o.b() ? _o.b()->getId() : -1)
    << ")";
  return _s;
}






/**********************************************************************************
 * ################################################################################
 * Boxed Objects Container supporting data.
 */



struct AxisOrderElm {

  AxisOrderElm() : aYAxisSearch(false), aMaxSearch(false), aBoxObjPtr(NULL) {
  };
  AxisOrderElm(bool _YAxis, bool _Max, Box *_ObjPtr) :
    aYAxisSearch(_YAxis),
    aMaxSearch(_Max),
    aBoxObjPtr(_ObjPtr) {
    };

  inline bool isXSearch() const { return aYAxisSearch; };
  inline bool isMaxSearch() const { return aMaxSearch; };

  inline int getId() const { return aBoxObjPtr ? aBoxObjPtr->getId() : -1; };

  inline bool operator < (AxisOrderElm const &_E) const {
#ifndef NDEBUG
    assert(aBoxObjPtr != NULL && _E.aBoxObjPtr != NULL);
#endif
    int i = (aYAxisSearch << 1) + aMaxSearch;
    int j = (_E.aYAxisSearch << 1) + _E.aMaxSearch;
    double &v1 = aBoxObjPtr->aBorder[i];
    double &v2 = _E.aBoxObjPtr->aBorder[j];
    return v1!=v2 ? v1<v2 : getId()<_E.getId(); // no 0 size particles !
  };

  inline bool operator > (AxisOrderElm const &_E) const {
    return ! operator < (_E);
  };


  bool aYAxisSearch : 1;
  bool aMaxSearch : 1;
  Box *aBoxObjPtr;


  friend ostream & operator << (ostream &_, AxisOrderElm const &_o);
};


ostream & operator << (ostream &_s, AxisOrderElm const &_o) {
  return _s
    << "{" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? "Max" : "Min")
    << ", Id:" << _o.getId()
    << "}";
}



/** Comparator for std::set to order axis.
*/
/*
struct AxisCmp {
  inline bool operator ()(AxisOrderElm const &_A, AxisOrderElm const &_B) const {
    return _A < _B;
  };
};
*/
// typedef set<AxisOrderElm, AxisCmp> T_AxisOrdered;
typedef set<AxisOrderElm> T_AxisOrdered;







/**********************************************************************************
 * ################################################################################
 * Boxed Objects Container and supporting data.
 */


typedef deque<BoxedObj> T_BoxPartCntr;
typedef deque<BoxedEdge> T_BoxEdgeCntr;

/**
*/
class BoxedObjCntr {

  public:

    BoxedObjCntr() {};

    inline bool isInit() const { return aInitialized; };

    void boxObjs(TwoDScene const &);
    void updateScene(TwoDScene const &);
    bool findAxisIntersect(T_AxisOrdered &, T_IntersSet &);
    bool findIntersect(T_IntersSet &);

    void update(TwoDScene const &, VectorXs const &, BoxedObj &);

    inline void addToOrderedAxis(BoxedObj &_O) {
      aXAxis.insert(AxisOrderElm(false, false, &_O));
      aXAxis.insert(AxisOrderElm(false, true, &_O));
      aYAxis.insert(AxisOrderElm(true, false, &_O));
      aYAxis.insert(AxisOrderElm(true, true, &_O));
    };


    //friend ostream & operator << (ostream &_, BoxedObj const &_o);


  private:

    inline void eraseFromOrderedAxis(BoxedObj &_O) {
      aXAxis.erase(AxisOrderElm(false, false, &_O));
      aXAxis.erase(AxisOrderElm(false, true, &_O));
      aYAxis.erase(AxisOrderElm(true, false, &_O));
      aYAxis.erase(AxisOrderElm(true, true, &_O));
    };


  private:

    bool aInitialized;
    double aAvgSize;
    T_BoxPartCntr aParticleBoxes; 
    T_BoxEdgeCntr aEdgeBoxes; 
    T_AxisOrdered aXAxis;
    T_AxisOrdered aYAxis;
};



/**
 * Boxes scene objects.
 */
void BoxedObjCntr::boxObjs(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  aAvgSize = 0.0;

  aParticleBoxes.resize(ParticlesNum);
  T_BoxPartCntr::iterator pIter = aParticleBoxes.begin();
  for (int i=0; i < ParticlesNum; ++i) {
    *pIter = BoxedObj(i, X.segment<2>(i<<1), _Scene.getRadius(i));
    addToOrderedAxis(*pIter);
    pIter++;
  };

  aAvgSize /= ParticlesNum;

  aInitialized = true;


#ifndef NDEBUG

  cout << __FUNCTION__
    << ", Particles=" << ParticlesNum
    << ", PBoxSize=" << aParticleBoxes.size()
    << ", AvgSize=" << aAvgSize
    << ", Boxes10:";
  // crashes! copy_n(aParticleBoxes.begin(), g_MaxCoutNum, ostream_iterator<BoxedObj>(cout, ", "));
  T_BoxPartCntr::const_iterator I=aParticleBoxes.begin();
  for (int i=0; i<g_MaxCoutNum && I!=aParticleBoxes.end(); ++i, ++I) {
    cout << *I << ", ";
  };

  cout << "\nXAxis:";
  T_AxisOrdered::const_iterator IA=aXAxis.begin(); 
  for(int i=0; i<g_MaxCoutNum && IA!=aXAxis.end(); ++IA) {
    cout << *IA << ", ";
  }
  cout << "\nYAxis:";
  IA=aYAxis.begin(); 
  for(int i=0; i<g_MaxCoutNum && IA!=aYAxis.end(); ++IA) {
    cout << *IA << ", ";
  }
  cout << endl;

#endif
}


void BoxedObjCntr::update(TwoDScene const &_Scene, VectorXs const &_X, BoxedObj &_O) {
  Vector2s X = _X.segment<2>(_O.getId()<<1);
  if (_O.changedX(X)) {

#ifndef NDEBUG
#ifdef MY_DEBUG
    cout << __FUNCTION__
      << " Detected particle to update:" << _O
      << ", with:" << X.transpose()
      << endl;
#endif
#endif

    eraseFromOrderedAxis(_O);
    _O.update(X, _Scene.getRadius(_O.getId()));
    addToOrderedAxis(_O);
  };
};



/**
 * Updates data from scene changes.
 */
void BoxedObjCntr::updateScene(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  T_BoxPartCntr::iterator pIter = aParticleBoxes.begin();
  for (; pIter != aParticleBoxes.end(); ++pIter) {
    update(_Scene, X, *pIter);
  };
}



bool BoxedObjCntr::findAxisIntersect(T_AxisOrdered &_Axis, T_IntersSet &_Intersects)
{
  bool bFound = false;
  set<int> Active;

  struct IntersectsInserter {
    IntersectsInserter(int _Id, T_IntersSet &_Collection) : aId(_Id), aCollection(_Collection) {
    };
    void operator ()(int _Id) {
      aCollection.insert(aId < _Id ? pair<int,int>(aId, _Id) : pair<int,int>(_Id, aId));
    };

    int aId;
    T_IntersSet & aCollection;
  };

  for (T_AxisOrdered::const_iterator I=_Axis.begin(); I != _Axis.end(); ++I) {
    int i = (*I).getId();

    if ((*I).isMaxSearch()) {
      Active.erase(i);
      continue;
    };

    if (Active.empty()) {
      Active.insert(i);
      continue;
    };

    IntersectsInserter Inserter(i, _Intersects); 
    for_each(Active.begin(), Active.end(), Inserter);

    Active.insert(i);

    bFound = true;
  };

  return bFound;
}



bool BoxedObjCntr::findIntersect(T_IntersSet &_Intersects)
{
  size_t OrigSize = _Intersects.size();
  T_IntersSet CandidatesX, CandidatesY;
  if (findAxisIntersect(aXAxis, CandidatesX)) {
    // Improvement: only search found candidates
    if (findAxisIntersect(aYAxis, CandidatesY)) {
      set_intersection(
          CandidatesX.begin(), CandidatesX.end(),
          CandidatesY.begin(), CandidatesY.end(),
          inserter(_Intersects, _Intersects.begin()));
    };
  };

  return _Intersects.size() > OrigSize;
}



static BoxedObjCntr g_PPCntr;




/**********************************************************************************
 * ################################################################################
 * Region.
 */
/*
   class Region : public Box {

   public:

   static void ConstructRegions(T_BoxPartCntr &_BoxedObjs);

   Region() {
   };


   public:
   deque<BoxedObj> aMembers;

   static deque<Region> cRegions;
   };
   */



/**********************************************************************************
 * ################################################################################
 * Local module functions.
 */



/*
   inline void boxP(Vector2s const &_X, int _R, Box &_B)
   {
   _B.first.first(_X[0]-_R);
   _B.first.second(_X[0]+_R);
   _B.second.first(_X[1]-_R);
   _B.second.second(_X[1]+_R);
   }
   */



static void estimateEuler(TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
  g_PPCntr.updateScene(_Scene);
  g_PPCntr.findIntersect(_PP);
}



static void estimateLagr(TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
}



/**
 * Detects potential collision pairs as all in.
 * Mostly for tests
 */
static void allIn(TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
  int ParticlesNum = _Scene.getNumParticles();
  int EdgesNum = _Scene.getNumEdges();
  int HalfPlanesNum = _Scene.getNumHalfplanes();

  combination(ParticlesNum, _PP);
  combination2(ParticlesNum, EdgesNum, _PE);
  combination2(ParticlesNum, HalfPlanesNum, _PH);
}



/**
 * Initializes the scene for E optimization.
 */
static void initEuler(TwoDScene const &_Scene)
{

}



/**
 * Initializes the scene for L optimization.
 */
static void initLagr(TwoDScene const &_Scene)
{
}




/**
 * Evaluates the scene to pre initialize data structures and select the best algorithm.
 */
static AlgoEnum evaluateScene(TwoDScene const &_Scene)
{
  if (!g_PPCntr.isInit()) {
    g_PPCntr.boxObjs(_Scene);
    initEuler(_Scene);
  };

  return EulerAlgo;
}




/**********************************************************************************
 * ################################################################################
 * ostream operators
 */



#ifndef NDEBUG

static ostream & operator << (ostream &_s, T_IntersIds const &_o) {
  return _s << "(" << _o.first << "," << _o.second << ")";
};


static ostream & operator << (ostream &_s, T_IntersSet const &_o) {
  // copy_n(_o.begin(), 10, ostream_iterator<T_IntersIds>(_s, ", "));
  T_IntersSet::const_iterator I = _o.begin();
  for (int i=0; i < g_MaxCoutNum && I != _o.end(); ++i, ++I) {
    _s << *I;
  };

  return _s;
};

#endif






/**********************************************************************************
 * ################################################################################
 * Global API implementation
 */


// Given particle positions, computes lists of *potentially* overlapping object
// pairs. How exactly to do this is up to you.
// Inputs: 
//   scene:  The scene object. Get edge information, radii, etc. from here. If 
//           for some reason you'd also like to use particle velocities in your
//           algorithm, you can get them from here too.
//   x:      The positions of the particle.
// Outputs:
//   pppairs: A list of (particle index, particle index) pairs of potentially
//            overlapping particles. IMPORTANT: Each pair should only appear
//            in the list at most once. (1, 2) and (2, 1) count as the same 
//            pair.
//   pepairs: A list of (particle index, edge index) pairs of potential
//            particle-edge overlaps.
//   phpairs: A list of (particle index, halfplane index) pairs of potential
//            particle-halfplane overlaps.
void ContestDetector::findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
  switch(evaluateScene(scene)) {
    case EulerAlgo:
      estimateEuler(scene, pppairs, pepairs, phpairs);
      break;

    case LagrAlgo:
      estimateLagr(scene, pppairs, pepairs, phpairs);
      break;

    case AllInAlgo:
    default:
      allIn(scene, pppairs, pepairs, phpairs);
  };

#ifndef NDEBUG

  cout << __FUNCTION__ << ": PP:" << pppairs << endl;

#endif
}
