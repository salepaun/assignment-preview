/**
 */

#include <iostream>
#include <set>
#include <utility>
#include <deque>
#include <algorithm>
#include <iterator>
#include <array>
#include <map>
#include <vector>

#ifndef NDEBUG
#include <assert.h>
#endif

#include "TwoDScene.h"
#include "ContestDetector.h"
#include "Ops.h"


using namespace std;


/**
 * ToDo Improvements:
 * 1. Search Y axis for only found candidates at X axis.
 * 2. Regions
 * 3. Auto method detection
 * 4. L. optimization
 * 5. Possibly only rechecking regions with moved (changed) objects at X and Y axises.
 * 6. Undetermined: Report touching objects (?)
 * 7. Optimize checking if particles belong to an edge.
 * 8. Optimize updates with checking change (X,Y) and updating only proper axis (X,Y).
 */


// For now: 0,1,2
#define MY_DEBUG 0


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


enum BoxType {
  UnknownBox = 0,
  VrtxBox = 1,
  EdgeBox = 2,
  HalpPlaneBox = 4,
  RegionBox = 8
};


/**********************************************************************************
 * ################################################################################
 * Local module variables.
 */


int g_MaxCoutNum = 10;



/**********************************************************************************
 * ################################################################################
 * Box.
 */



class Box {

  public:

    typedef pair<int, BoxType> T_Key;
    typedef map<T_Key, Box*> T_Parent;
    typedef pair<T_Parent::key_type, T_Parent::mapped_type> T_ParentElm;

    static char const * typeAsStr(BoxType _Type) {
      switch (_Type) {
        case VrtxBox: return "V";
        case EdgeBox: return "E";
        case HalpPlaneBox: return "H";
        default: return "x";
      };
    };

    Box() : aId(-1), aType(UnknownBox), aChanged(false), aInitialized(false) {};
    Box(int _Id, BoxType _Type) : aId(_Id), aType(_Type), aChanged(false), aInitialized(false) {};

    inline int getId() const { return aId; };
    inline BoxType getType() const { return aType; };
    inline T_Key getKey() const { return T_Key(aId, aType); };

    inline bool isInitialized() const { return aInitialized; };
    inline void setInitialized() { aInitialized = true; };

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

    inline void registerParent(Box &_Parent) {
      aParents.insert(T_ParentElm(_Parent.getKey(), &_Parent));
    };
    inline void deregisterParent(Box &_Parent) {
      aParents.erase(_Parent.getKey());
    };
    inline bool hasParent(Box const &_Parent) const {
      return aParents.find(_Parent.getKey()) != aParents.end();
    };

    inline bool hasChanged() const {
      return aChanged;
    };
    inline void resetChange() {
      aChanged = false;
    };

    void changed() {
      aChanged = true;
      for_each(aParents.begin(), aParents.end(), 
          [](T_ParentElm const &_p) { _p.second->changed(); });
    };


    friend ostream & operator << (ostream &_, Box const &_o);

  public:

    int aId;
    BoxType aType;
    bool aChanged;
    bool aInitialized;
    T_Parent aParents;
    array<double, 4> aBorder;
};




ostream & operator << (ostream &_s, Box::T_Key const &_o) {
  return _s << _o.first << ":" << Box::typeAsStr(_o.second);
}


ostream & operator << (ostream &_s, Box const &_o) {
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
  _s << "#{" << _o.aId << ":" << Box::typeAsStr(_o.aType)
    << " (mnX:" << _o.aBorder[0]
    << ", mxX:" << _o.aBorder[1]
    << ", mnY:" << _o.aBorder[2]
    << ", mxY:" << _o.aBorder[3]
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
    BoxedObj(int _Id) : Box(_Id, VrtxBox) {};
    BoxedObj(int _Id, bool _Fixed, Vector2s const &_X, double const &_R) : 
      Box(_Id, VrtxBox), aFixed(_Fixed), aX(_X[0]), aY(_X[1]) {
      // !!! aAvgSize += _R;
      updateR(_R);
    };

    /** 
     * Returns true for fixed vertexes.
     *
     * Similar functionality for edges - no virtual methods in the box though -
     * high performance penalty for virtual objects - to improve if otherwise.
     */
    inline bool isFixed() const {
      return aFixed;
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
      changed();
    };


    friend ostream & operator << (ostream &_, BoxedObj const &_o);


  private:

    bool aFixed;
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
    BoxedEdge(int _Id) : Box(_Id, EdgeBox), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id, BoxedObj *_A, BoxedObj *_B, double const &_R) : 
      Box(_Id, EdgeBox), aVrtx(_A, _B) {
        updateR(_R);
      };


    inline bool initialized() const {
      return aVrtx.first && aVrtx.second;
    };

    inline void registerSelf() {
      if (initialized()) {
        aVrtx.first->registerParent(*this);
        aVrtx.second->registerParent(*this);
      };
    };

    /** 
     * Returns true for fixed vertexes.
     *
     * Similar functionality for edges - no virtual methods in the box though -
     * high performance penalty for virtual objects - to improve if otherwise.
     */
    inline bool isFixed() const {
      return initialized() && (a()->isFixed() && b()->isFixed());
    };

    inline bool isMember(int _Id) const {
      return initialized() && (a()->getId() == _Id || b()->getId() == _Id);
    };

    inline bool changedX() const {
      // should update by it's children
      // initialized() && (a()->hasChanged() || b()->hasChanged());
      return hasChanged();
    };

    inline void updateR(double const &_R) {
      if (initialized()) {
        aBorder[0] = (a()->xmin() < b()->xmin() ? a()->xmin() : b()->xmin()) - _R;
        aBorder[1] = (a()->xmax() > b()->xmax() ? a()->xmax() : b()->xmax()) + _R;
        aBorder[2] = (a()->ymin() < b()->ymin() ? a()->ymin() : b()->ymin()) - _R;
        aBorder[3] = (a()->ymax() > b()->ymax() ? a()->ymax() : b()->ymax()) + _R;
      };
    };

    inline void update(double const &_R) {
      updateR(_R);
      changed();
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

    bool aFixed;
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

  AxisOrderElm() : aYAxisSearch(false), aMaxSearch(false), aBoxObjPtr(NULL)
  { };

  AxisOrderElm(bool _YAxis, bool _Max, Box *_ObjPtr) :
    aYAxisSearch(_YAxis),
    aMaxSearch(_Max),
    aBoxObjPtr(_ObjPtr)
  { };

  inline bool isXSearch() const { return aYAxisSearch; };
  inline bool isMaxSearch() const { return aMaxSearch; };

  inline int getId() const { return aBoxObjPtr ? aBoxObjPtr->getId() : -1; };
  inline BoxType getType() const { return aBoxObjPtr ? aBoxObjPtr->getType() : UnknownBox; };
  inline Box::T_Key getKey() const { return aBoxObjPtr ? aBoxObjPtr->getKey() : Box::T_Key(-1, UnknownBox); };

  inline bool operator < (AxisOrderElm const &_E) const {
#ifndef NDEBUG
    assert(aBoxObjPtr != NULL && _E.aBoxObjPtr != NULL);
#endif
    int i = (aYAxisSearch << 1) + aMaxSearch;
    int j = (_E.aYAxisSearch << 1) + _E.aMaxSearch;
    double &v1 = aBoxObjPtr->aBorder[i];
    double &v2 = _E.aBoxObjPtr->aBorder[j];
    return v1!=v2 ? v1<v2 : getId()!=_E.getId() ? getId()<_E.getId() : getType()<_E.getType();
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
    << ", Id:" << _o.getId() << ":" << Box::typeAsStr(_o.getType())
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
 * Scene region
 */


typedef map<int, BoxedObj *> T_BoxPartMap;
typedef pair<int, BoxedObj *> T_MapPair;


typedef vector<BoxedObj> T_BoxPartCntr;
typedef vector<BoxedEdge> T_BoxEdgeCntr;

typedef map<int, BoxedObj *> T_BoxPartMap;
typedef pair<int, BoxedObj *> T_MapPair;




/**
*/
class BoxedRegion : Box {

  public:

    BoxedRegion() : Box(), apBoxedParts(NULL), apBoxedEdges(NULL) {};
    BoxedRegion(int _Id, T_BoxPartCntr const *_BoxedParts, T_BoxEdgeCntr const *_BoxedEdges) : 
      Box(_Id, RegionBox), apBoxedParts(_BoxedParts), apBoxedEdges(_BoxedEdges) {
        setInitialized();
      };

    bool findAxisIntersect(T_AxisOrdered &, T_IntersSet &, T_IntersSet &, T_IntersSet &);
    bool findIntersect(T_IntersSet &, T_IntersSet &, T_IntersSet &);

    inline void add(Box &_B) {
      addToOrderedAxis(_B);
    };
    inline void erase(Box &_B) {
      eraseFromOrderedAxis(_B);
    };


    friend ostream & operator << (ostream &_, BoxedRegion const &_o);


  private:

    inline void addToOrderedAxis(Box &_O) {
      aXAxis.insert(AxisOrderElm(false, false, &_O));
      aXAxis.insert(AxisOrderElm(false, true, &_O));
      aYAxis.insert(AxisOrderElm(true, false, &_O));
      aYAxis.insert(AxisOrderElm(true, true, &_O));
    };

    inline void eraseFromOrderedAxis(Box &_O) {
      aXAxis.erase(AxisOrderElm(false, false, &_O));
      aXAxis.erase(AxisOrderElm(false, true, &_O));
      aYAxis.erase(AxisOrderElm(true, false, &_O));
      aYAxis.erase(AxisOrderElm(true, true, &_O));
    };


  private:

    T_AxisOrdered aXAxis;
    T_AxisOrdered aYAxis;

    T_BoxPartCntr const *apBoxedParts;
    T_BoxEdgeCntr const *apBoxedEdges;
};




ostream & operator << (ostream &_s, BoxedRegion const &_o) {
  _s << Box(_o);
  dumpContainer<>(g_MaxCoutNum, _s, NULL, "\n ****  XAxis:",
      _o.aXAxis.size(), _o.aXAxis.begin(), _o.aXAxis.end());

  dumpContainer<>(g_MaxCoutNum, _s, NULL, "\n ****  YAxis:",
      _o.aYAxis.size(), _o.aYAxis.begin(), _o.aYAxis.end());

  return _s;
}




bool BoxedRegion::findAxisIntersect(
    T_AxisOrdered &_Axis,
    T_IntersSet &_PP,
    T_IntersSet &_PE,
    T_IntersSet &_PH)
{
  bool bFound = false;
  set<Box::T_Key> Active;

  struct IntersectsInserter {
    IntersectsInserter(Box::T_Key _Key, 
        T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH,
        T_BoxEdgeCntr const *_pEdgeCntr) : 
      aKey(_Key), 
      aPP(_PP),
      aPE(_PE),
      aPH(_PH),
      apEdgeCntr(_pEdgeCntr) {
    };

    void operator ()(Box::T_Key _Key) {
      T_IntersSet *intersects = NULL;
      int Id1, Id2;
      switch((aKey.second | _Key.second)) {
        case VrtxBox: 
          {
            if (aKey.first < _Key.first) {
              Id1 = aKey.first; Id2 = _Key.first;
            } else {
              Id1 = _Key.first; Id2 = aKey.first;
            };
            intersects = &aPP;
          }; break;
        case EdgeBox:
        case 3:
          {
            if (aKey.second == EdgeBox) {
              Id1 = _Key.first; Id2 = aKey.first;
            } else {
              Id1 = aKey.first; Id2 = _Key.first;
            }
            intersects = apEdgeCntr ? apEdgeCntr->at(Id2).isMember(Id1) ?  NULL : &aPE : &aPE;
          }; break;
        case HalpPlaneBox:
        case 5:
          {
            if (aKey.second == HalpPlaneBox) {
              Id1 = _Key.first; Id2 = aKey.first;
            } else {
              Id1 = aKey.first; Id2 = _Key.first;
            }
            // intersects = aEdgeCntr[Id2].isMember(Id1) ?  NULL : &aPE;
            intersects = &aPH; break;
          };
      };

#ifndef NDEBUG
#ifdef MY_DEBUG
#if MY_DEBUG > 2
      cout << __FUNCTION__
        << ", Akey:" << aKey
        << ", Bkey:" << _Key
        << ", Atype|Btype:" << (aKey.second|_Key.second)
        << ", Intersecs:" << intersects
        << endl;
#endif
#endif
#endif
      if(intersects) {
        intersects->insert(pair<int,int>(Id1, Id2));
      };
    };

    Box::T_Key aKey;
    T_IntersSet &aPP, &aPE, &aPH;
    T_BoxEdgeCntr const *apEdgeCntr;
  };

  for (T_AxisOrdered::const_iterator I=_Axis.begin(); I != _Axis.end(); ++I) {
    Box::T_Key key = (*I).getKey();


#ifndef NDEBUG
#ifdef MY_DEBUG
#if MY_DEBUG > 2
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Active:",
        Active.size(), Active.begin(), Active.end()) << endl;
#endif
#endif
#endif

    if ((*I).isMaxSearch()) {
      Active.erase(key);
      continue;
    };

    if (Active.empty()) {
      Active.insert(key);
      continue;
    };

    IntersectsInserter Inserter(key, _PP, _PE, _PH, apBoxedEdges);
    for_each(Active.begin(), Active.end(), Inserter);

    Active.insert(key);

    bFound = true;
  };

  return bFound;
}



bool BoxedRegion::findIntersect(
    T_IntersSet &_PP,
    T_IntersSet &_PE,
    T_IntersSet &_PH)
{
  size_t PPSize = _PP.size();
  size_t PESize = _PE.size();
  size_t PHSize = _PH.size();

  T_IntersSet XInterPP, XInterPE, XInterPH, YInterPP, YInterPE, YInterPH;
  if (findAxisIntersect(aXAxis, XInterPP, XInterPE, XInterPH)) {
    // Improvement: only search found candidates !!!
    if (findAxisIntersect(aYAxis, YInterPP, YInterPE, YInterPH)) {
      if (XInterPP.size() && YInterPP.size()){
        set_intersection(
            XInterPP.begin(), XInterPP.end(),
            YInterPP.begin(), YInterPP.end(),
            inserter(_PP, _PP.begin()));
      };

      if (XInterPE.size() && YInterPE.size()){
        set_intersection(
            XInterPE.begin(), XInterPE.end(),
            YInterPE.begin(), YInterPE.end(),
            inserter(_PE, _PE.begin()));
      };

      if (XInterPH.size() && YInterPH.size()){
        set_intersection(
            XInterPH.begin(), XInterPH.end(),
            YInterPH.begin(), YInterPH.end(),
            inserter(_PH, _PH.begin()));
      };
    };
  };

#ifndef NDEBUG
#ifdef MY_DEBUG
  cout << __FUNCTION__
    << " PPx:" << XInterPP.size()
    << " PEx:" << XInterPE.size()
    << " PHx:" << XInterPH.size()
    << " PPy:" << YInterPP.size()
    << " PEy:" << YInterPE.size()
    << " PHy:" << YInterPH.size()
    << endl;
#endif
#endif

  return _PP.size() > PPSize || _PE.size() > PESize || _PH.size() > PHSize;
}






/**********************************************************************************
 * ################################################################################
 * Boxed Objects Container and supporting data.
 */


/**
*/
class BoxedScene {

  public:

    BoxedScene() : aRegion(0, &aBoxParts, &aBoxEdges) {};

    inline bool isInit() const { return aInitialized; };

    inline void init(TwoDScene const &_Scene) {
      boxObjs(_Scene);
    };
    void update(TwoDScene const &);

    inline bool findIntersect(T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) {
      return aRegion.findIntersect(_PP, _PE, _PH);
    };

    inline void addToRegion(Box &_Box) {
      aRegion.add(_Box);
    };
    inline void eraseFromRegion(Box &_Box) {
      aRegion.erase(_Box);
    };

    //friend ostream & operator << (ostream &_, BoxedObj const &_o);


  private:

    void boxVrtxs(TwoDScene const &);
    void boxEdges(TwoDScene const &);
    void boxObjs(TwoDScene const &);

    void updateVrtx(TwoDScene const &, VectorXs const &, BoxedObj &);
    void updateEdge(TwoDScene const &, BoxedEdge &);


  private:

    bool aInitialized;
    double aAvgSize;

    BoxedRegion aRegion;

    T_BoxPartCntr aBoxParts; 
    T_BoxEdgeCntr aBoxEdges; 

    T_BoxPartMap aParticleMap;
};



/**
 * Boxes scene objects.
 */
void BoxedScene::boxVrtxs(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  aBoxParts.resize(ParticlesNum);
  T_BoxPartCntr::iterator pIter = aBoxParts.begin();
  for (int i=0; i < ParticlesNum; ++i) {
    *pIter = BoxedObj(i, _Scene.isFixed(i), X.segment<2>(i<<1), _Scene.getRadius(i));
    addToRegion(*pIter);
    aParticleMap.insert(T_MapPair(i, &(*pIter)));
    pIter++;
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxEdges(TwoDScene const &_Scene)
{
  int EdgesNum = _Scene.getNumEdges();

  aBoxEdges.resize(EdgesNum);
  T_BoxEdgeCntr::iterator pIter = aBoxEdges.begin();
  for (int i=0; i < EdgesNum; ++i) {
    double R = _Scene.getEdgeRadii()[i];
    T_IntersIds Vrtxs = _Scene.getEdge(i);
    *pIter = BoxedEdge(i, aParticleMap[Vrtxs.first], aParticleMap[Vrtxs.second], R);
    addToRegion(*pIter);
    pIter++;
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxObjs(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  aAvgSize = 0.0;

  boxVrtxs(_Scene);
  // Edges have to go after Vertexes!
  boxEdges(_Scene);

  aInitialized = true;


#ifndef NDEBUG

  cout << __FUNCTION__
    << ", Particles=" << ParticlesNum
    << ", PBoxSize=" << aBoxParts.size()
    << ", AvgSize=" << aAvgSize;

  // crashes! copy_n(aBoxParts.begin(), g_MaxCoutNum, ostream_iterator<BoxedObj>(cout, ", "));
  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** VrtxBox:",
      aBoxParts.size(), aBoxParts.begin(), aBoxParts.end());

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** EdgeBox:",
      aBoxEdges.size(), aBoxEdges.begin(), aBoxEdges.end());

  cout << aRegion
    << endl;

#endif
}


void BoxedScene::updateVrtx(TwoDScene const &_Scene, VectorXs const &_X, BoxedObj &_O) {
  Vector2s X = _X.segment<2>(_O.getId()<<1);
  if (_O.changedX(X)) {

#ifndef NDEBUG
#ifdef MY_DEBUG 
#if MY_DEBUG > 1
    cout << __FUNCTION__
      << " Detected particle to update:" << _O
      << ", with:" << X.transpose()
      << endl;
#endif
#endif
#endif

    eraseFromRegion(_O);
    _O.update(X, _Scene.getRadius(_O.getId()));
    addToRegion(_O);
  };
};


void BoxedScene::updateEdge(TwoDScene const &_Scene, BoxedEdge &_O) {
  if (_O.changedX()) {

#ifndef NDEBUG
#ifdef MY_DEBUG
#if MY_DEBUG > 1
    cout << __FUNCTION__
      << " Detected edge to update:" << _O
      << endl;
#endif
#endif
#endif

    eraseFromRegion(_O);
    _O.update(_Scene.getEdgeRadii()[_O.getId()]);
    addToRegion(_O);
  };
};



/**
 * Updates data from scene changes.
 */
void BoxedScene::update(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  T_BoxPartCntr::iterator pIter = aBoxParts.begin();
  for (; pIter != aBoxParts.end(); ++pIter) {
    if (!(*pIter).isFixed()) {
      updateVrtx(_Scene, X, *pIter);
    };
  };

  // Edges have to go after Vertexes!
  T_BoxEdgeCntr::iterator pIterE = aBoxEdges.begin();
  for (; pIterE != aBoxEdges.end(); ++pIterE) {
    if (!(*pIter).isFixed()) {
      updateEdge(_Scene, *pIterE);
    };
  };

  for_each(aBoxParts.begin(), aBoxParts.end(), [](BoxedObj &_o) { _o.resetChange(); });
  for_each(aBoxEdges.begin(), aBoxEdges.end(), [](BoxedEdge &_o) { _o.resetChange(); });

#ifndef NDEBUG
#ifdef MY_DEBUG
#if MY_DEBUG > 1

  cout << aRegion << endl;

#endif
#endif
#endif
}



static BoxedScene g_Scene;




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



static void estimateEuler(
    bool _First,
    TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
  if (!_First) {
    g_Scene.update(_Scene);
  };
  g_Scene.findIntersect(_PP, _PE, _PH);
}



static void estimateLagr(
    bool _First,
    TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
}



/**
 * Detects potential collision pairs as all in.
 * Mostly for tests
 */
static void allIn(
    bool _First,
    TwoDScene const &_Scene,
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
static AlgoEnum evaluateScene(TwoDScene const &_Scene, bool &_First)
{
  _First = !g_Scene.isInit();
  if (_First) {
    g_Scene.init(_Scene);
    initEuler(_Scene);
  };

  return AllInAlgo;
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
  return dumpContainer<>(g_MaxCoutNum, _s, NULL, NULL, _o.size(), _o.begin(), _o.end());
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
  bool bFirst = false;

  switch(evaluateScene(scene, bFirst)) {
    case EulerAlgo:
      estimateEuler(bFirst, scene, pppairs, pepairs, phpairs);
      break;

    case LagrAlgo:
      estimateLagr(bFirst, scene, pppairs, pepairs, phpairs);
      break;

    case AllInAlgo:
    default:
      allIn(bFirst, scene, pppairs, pepairs, phpairs);
  };

#ifndef NDEBUG
#ifdef MY_DEBUG

  if (pppairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** PP:",
        pppairs.size(), pppairs.begin(), pppairs.end()) << endl;
  if (pepairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** PE:",
        pepairs.size(), pepairs.begin(), pepairs.end()) << endl;
  if (phpairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** PH:",
        phpairs.size(), phpairs.begin(), phpairs.end()) << endl;

#endif
#endif
}
