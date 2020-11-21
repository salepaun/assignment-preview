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
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <ctime>

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
 * 9. To replace axises sets with deque/vector and to sort myself after update and before finding intersections
 * 10. Save intersections and use if no changes, discard on changes.
 * 11. Copy ordered axises to children after split instead of reinserting.
 * 12. Keep references to particles in ordered axises and operate only if needed (1 check for 4 set operations).
 */


// For now: 0,1,2,3,4,5
#define MY_DEBUG 0
// 0,1,2
#define MY_TIMING 1


#define D_TIME_SEC(a,b) cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< endl; a = clock();
#define D_TIME_SEC2(a,b,c) cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< endl; a = clock();
#define D_TIME_SEC3(a,b,c,d) cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< ": " << d \
<< endl; a = clock();
#define D_TIME_SEC4(a,b,c,d,e) cout << "= TIME:" << __FUNCTION__ \
  << ": " << (float)(clock()-a)/CLOCKS_PER_SEC << "[s]" \
<< ": " << b \
<< ": " << c \
<< ": " << d \
<< ": " << e \
<< endl; a = clock();


/**********************************************************************************
 * ################################################################################
 * Local module declarations.
 */


typedef pair<double, double> T_MinMax;
typedef pair<int, int> T_IntersIds;
typedef vector<T_IntersIds> T_IntersCntr;


enum AlgoEnum {
  AllInAlgo = 0,
  EulerAlgo,
  LagrAlgo
};


enum BoxType {
  E_UnknownBox = 0,
  E_VrtxBox = 1,
  E_EdgeBox = 2,
  E_HalfPBox = 4,
  E_RegionBox = 8
};


/**********************************************************************************
 * ################################################################################
 * Local module variables.
 */


int g_MaxCoutNum = 20;



/**********************************************************************************
 * ################################################################################
 * Box.
 */



ostream & operator << (ostream &_s, array<double,4> const &_o) {
  copy(_o.begin(), _o.end(), ostream_iterator<double>(_s, ","));
  return _s;
}


class Box
{

  public:

    typedef pair<int, BoxType> T_Key;
    typedef map<T_Key, Box*> T_Parent;
    typedef pair<T_Parent::key_type, T_Parent::mapped_type> T_ParentElm;
    typedef array<double, 4> T_Boundry;
    typedef vector<Box const *> T_BoxCntr;


    /**
     * Hasher.
     * Needed for unordered_set.
     */
    struct KeyHash {
      inline size_t operator ()(T_Key const &_K) const {
        return hash<int>{}(_K.first) ^ hash<int>{}(_K.second);
      };
    };

    static char const * typeAsStr(BoxType _Type) {
      switch (_Type) {
        case E_VrtxBox: return "V";
        case E_EdgeBox: return "E";
        case E_HalfPBox: return "H";
        case E_RegionBox: return "R";
        default: return "x";
      };
    };

    Box() : 
      aId(-1), aType(E_UnknownBox),
      aChanged(false), aInitialized(false), aFixed(false) {};
    Box(int _Id, BoxType _Type, bool _Fixed=false) :
      aId(_Id), aType(_Type),
      aChanged(false), aInitialized(false), aFixed(_Fixed) {};

    inline int getId() const { return aId; };
    inline BoxType getType() const { return aType; };
    inline T_Key getKey() const { return T_Key(aId, aType); };

    inline bool isInitialized() const { return aInitialized; };
    inline void setInitialized() {
      aInitialized = true; };

    inline bool isFixed() const { return aFixed; };
    inline void setFixed(bool _Fixed=true) { aFixed = _Fixed; };

    inline bool hasXInf() const { return (std::isinf(aBorder[0]) || std::isinf(aBorder[1])); };
    inline bool hasYInf() const { return (std::isinf(aBorder[2]) || std::isinf(aBorder[3])); };
    inline bool hasInf() const { return (hasXInf() || hasYInf()); };

    inline double const & xmin() const { return aBorder[0]; };
    inline double const & xmax() const { return aBorder[1]; };
    inline double const & ymin() const { return aBorder[2]; };
    inline double const & ymax() const { return aBorder[3]; };

    inline void setXMin(double const &_Val) { aBorder[0] = _Val; };
    inline void setXMax(double const &_Val) { aBorder[1] = _Val; };
    inline void setYMin(double const &_Val) { aBorder[2] = _Val; };
    inline void setYMax(double const &_Val) { aBorder[3] = _Val; };

    inline double xsize() const { return aBorder[1] - aBorder[0]; };
    inline double ysize() const { return aBorder[3] - aBorder[2]; };

    inline void shrinkXBy(double _XFactor) { 
      aBorder[1] = aBorder[0] + xsize() / abs(_XFactor);
    };
    inline void shrinkYBy(double _YFactor) {
      aBorder[3] = aBorder[2] + ysize() / abs(_YFactor);
    };
    inline void shrinkBy(double _Factor) {
      shrinkXBy(_Factor); shrinkYBy(_Factor);
    };

    inline bool safeShiftX(bool _Right=true) {
      if (hasXInf())
        return false;
      double XSize = xsize();
      if (_Right) {
        aBorder[0] = aBorder[1];
        aBorder[1] += XSize;
      } else {
        aBorder[1] = aBorder[0];
        aBorder[0] -= XSize;
      };
      return true;
    };

    inline bool safeShiftY(bool _Right=true) {
      if (hasYInf())
        return false;
      double YSize = ysize();
      if (_Right) {
        aBorder[2] = aBorder[3];
        aBorder[3] += YSize;
      } else {
        aBorder[3] = aBorder[2];
        aBorder[2] -= YSize;
      };
      return true;
    };

    inline bool intersects(Box const &_B) const {
      return intersects(_B.getBoundry());
    };
    inline bool intersects(T_Boundry const &_B) const {
      bool R = intersectsAxis(_B[0], _B[1], aBorder[0], aBorder[1])
        && intersectsAxis(_B[2], _B[3], aBorder[2], aBorder[3]);
      return R;
    };
    inline bool intersectsX(Box const &_B) const {
      return intersectsAxis(_B.aBorder[0], _B.aBorder[1], aBorder[0], aBorder[1]);
    };
    inline bool intersectsY(Box const &_B) const {
      return intersectsAxis(_B.aBorder[2], _B.aBorder[3], aBorder[2], aBorder[3]);
    };
    inline bool intersectsAxis(double const &_An, double const &_Ax,
        double const &_Bn, double const &_Bx) const {
      return !(_An > _Bx || _Ax < _Bn);
    };
    inline bool intersectsXOrdered(Box const &_B) const {
      return aBorder[1] > _B.aBorder[0];
    };
    inline bool intersectsYOrdered(Box const &_B) const {
      return aBorder[3] > _B.aBorder[2];
    };
    inline bool intersectsAxisOld(double const &_An, double const &_Ax,
        double const &_Bn, double const &_Bx) const {
      return (_An <= _Bx && _An >= _Bn) || (_Bn <= _Ax && _Bn >= _An);
    };

    inline void expandNoInf(Box const &_B) {
      expandNoInf(_B.getBoundry());
    };
    inline void expandNoInf(T_Boundry const &_B) {
      if (!std::isinf(_B[0]) && (_B[0] < aBorder[0] || std::isinf(aBorder[0]))) aBorder[0] = _B[0];
      if (!std::isinf(_B[1]) && (_B[1] > aBorder[1] || std::isinf(aBorder[1]))) aBorder[1] = _B[1];
      if (!std::isinf(_B[2]) && (_B[2] < aBorder[2] || std::isinf(aBorder[2]))) aBorder[2] = _B[2];
      if (!std::isinf(_B[3]) && (_B[3] > aBorder[3] || std::isinf(aBorder[3]))) aBorder[3] = _B[3];
    };
    inline void expand(Box const &_B) {
      expand(_B.getBoundry());
    };
    inline void expand(T_Boundry const &_B) {
      if (_B[0] < aBorder[0]) aBorder[0] = _B[0];
      if (_B[1] > aBorder[1]) aBorder[1] = _B[1];
      if (_B[2] < aBorder[2]) aBorder[2] = _B[2];
      if (_B[3] > aBorder[3]) aBorder[3] = _B[3];
    };

    inline T_Boundry const & getBoundry() const { return aBorder; };
    inline void setBoundry(Box const &_B) {
      setBoundry(_B.getBoundry());
    };
    inline void setBoundry(double const (&_B)[4]) {
      // copy(aBorder.begin(), aBorder.end(), std::begin(_B));
      memcpy(aBorder.data(), _B, sizeof(_B));
    };
    inline void setBoundry(T_Boundry const &_B) {
      aBorder = _B;
    };
    inline void setBoundry(double const &_Xmn, double const &_Xmx,
        double const &_Ymn, double const &_Ymx) {
      aBorder = {_Xmn, _Xmx, _Ymn, _Ymx};
    };
    inline void setMinMaxBoundry() {
      aBorder = {
        numeric_limits<double>::min(),
        numeric_limits<double>::max(),
        numeric_limits<double>::min(),
        numeric_limits<double>::max()};
    };
    inline void setInfBoundry() {
      aBorder = {
        -numeric_limits<double>::infinity(),
        numeric_limits<double>::infinity(),
        -numeric_limits<double>::infinity(),
        numeric_limits<double>::infinity()};
    };

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
      if (!hasParent(_Parent)) {
        aParents.insert(T_ParentElm(_Parent.getKey(), &_Parent));
        _Parent.changed();
      };
    };
    inline void deregisterParent(Box &_Parent) {
      if (hasParent(_Parent)) {
        aParents.erase(_Parent.getKey());
        _Parent.changed();
      };
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
          [](T_ParentElm const &_p) { if(_p.second) _p.second->changed(); });
    };


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, Box const &);


  public:

    int aId;
    BoxType aType;
    bool aChanged;
    bool aInitialized;
    bool aFixed;
    T_Parent aParents;
    T_Boundry aBorder;
};




ostream & operator << (ostream &_s, Box::T_Key const &_o) {
  return _s << _o.first << ":" << Box::typeAsStr(_o.second);
}


ostream & Box::toStr(ostream &_s) const {
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
  return _s << "#{" << aId << ":" << Box::typeAsStr(aType)
    << " Init:" << aInitialized
    << " PNum:" << aParents.size()
    << " Chng:" << aChanged
    << " (mnX:" << aBorder[0]
    << ", mxX:" << aBorder[1]
    << ", mnY:" << aBorder[2]
    << ", mxY:" << aBorder[3]
    << ")";
}


ostream & operator << (ostream &_s, Box const &_o) {
  return _o.toStr(_s);
}


/**********************************************************************************
 * ################################################################################
 * Boxed Object.
 */

/**
*/
class BoxedVrtx : public Box
{

  public:

    BoxedVrtx() : Box() {};
    BoxedVrtx(int _Id) : Box(_Id, E_VrtxBox) {};
    BoxedVrtx(int _Id, bool _Fixed, Vector2s const &_X, double const &_R) : 
      Box(_Id, E_VrtxBox, _Fixed), aX(_X[0]), aY(_X[1]) {
        // !!! aAvgSize += _R;
        updateR(_R);
        changed();
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


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedVrtx const &);


  private:

    double aX;
    double aY;
};



ostream & BoxedVrtx::toStr(ostream &_s) const {
  return Box::toStr(_s)
    << " (" << aX << "," << aY
    << ")";
}


ostream & operator << (ostream &_s, BoxedVrtx const &_o) {
  return _o.toStr(_s);
};



/**********************************************************************************
 * ################################################################################
 * Boxed Edge.
 */

/**
*/
class BoxedEdge : public Box
{
  public:

    typedef pair<BoxedVrtx *, BoxedVrtx *> T_VrtxPair;

    BoxedEdge() : Box(), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id) : Box(_Id, E_EdgeBox, false), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id, BoxedVrtx *_A, BoxedVrtx *_B, double const &_R) : 
      Box(_Id, E_EdgeBox, false), aVrtx(_A, _B) {
        updateR(_R);
        changed();
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
     * Updates it's fixed flag.
     *
     * Similar functionality for edges - no virtual methods in the box though -
     * high performance penalty for virtual objects - to improve if otherwise.
     */
    inline void setEdgeFixed() {
      setFixed(initialized() && (a()->isFixed() && b()->isFixed()));
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


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedEdge const &);


  private:

    inline BoxedVrtx const * a() const {
      return aVrtx.first;
    };
    inline BoxedVrtx const * b() const {
      return aVrtx.second;
    };


  private:

    T_VrtxPair aVrtx;
};



ostream & BoxedEdge::toStr(ostream &_s) const {
  return Box::toStr(_s)
    << " (A:" << (a() ? a()->getId() : -1)
    << ",B:" << (b() ? b()->getId() : -1)
    << ")";
}


ostream & operator << (ostream &_s, BoxedEdge const &_o) {
  return _o.toStr(_s);
}




/**********************************************************************************
 * ################################################################################
 * Boxed Halfplane.
 */

/**
 * Boxed halfplane of the scene.
 */
class BoxedHalfP : public Box
{
  public:

    BoxedHalfP() : Box() {};
    BoxedHalfP(int _Id) : Box(_Id, E_HalfPBox, true) {};
    BoxedHalfP(int _Id, Vector2s const &_X, Vector2s const &_N) :
      Box(_Id, E_HalfPBox, true), aX(_X[0]), aY(_X[1]), aNx(_N[0]), aNy(_N[1]) {
        update();
        changed();
      };

    inline void update() {
      setInfBoundry();
      if (!aNx && aNy) aBorder[2 + (aNy>0)] = aY; 
      if (aNx && !aNy) aBorder[(aNx>0)] = aX; 
    };


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedVrtx const &);


  private:

    double aX, aY;
    double aNx, aNy;
};



ostream & BoxedHalfP::toStr(ostream &_s) const {
  return Box::toStr(_s)
    << " x(" << aX << "," << aY
    << "),n(" << aNx << "," << aNy << ")";
}


ostream & operator << (ostream &_s, BoxedHalfP const &_o) {
  return _o.toStr(_s);
};








/**********************************************************************************
 * ################################################################################
 * Boxed Objects Container supporting data.
 */



/**
 * Element for ordered axis.
 */
struct OrderedAxisElm
{
  OrderedAxisElm() : aKey(Box::T_Key(-1, E_UnknownBox)), 
  aYAxisSearch(false), aMaxSearch(false),
  aFixed(false),
  apBox(NULL)
  {
  };

  OrderedAxisElm(bool _YAxis, bool _Max, Box &_B) :
    aYAxisSearch(_YAxis),
    aMaxSearch(_Max),
    aFixed(false),
    apBox(&_B)
  { 
    aKey = _B.getKey();
    aFixed = _B.isFixed();
    aVal = getRawVal(*this, _B);
  };


  OrderedAxisElm const & initAxisPair(OrderedAxisElm const &_E) {
    aKey = _E.getKey();
    aFixed = _E.isFixed();
    apBox = _E.apBox;
    aMaxSearch = !_E.aMaxSearch;
    aYAxisSearch = _E.aYAxisSearch;
    aVal = getRawVal(*this, *apBox);
    return *this;
  };


  inline bool isXSearch() const { return aYAxisSearch; };
  inline bool isMaxSearch() const { return aMaxSearch; };

  inline int getId() const {
    return aKey.first;
  };
  inline BoxType getType() const {
    return aKey.second;
  };
  inline Box::T_Key getKey() const {
    return aKey;
  };
  inline bool isFixed() const {
    return aFixed;
  };
  inline double const & getVal() const {
    return aVal;
  };
  inline double const & getRawVal(OrderedAxisElm const &_E, Box &_B) const {
    int i = (_E.aYAxisSearch << 1) + _E.aMaxSearch;
    return _B.aBorder[i];
  };

  inline Box const * getBox() const { return apBox; };

  inline void changed() {
    if (apBox) apBox->changed();
  };
  inline void registerParent(Box &_Parent) {
    if (apBox) apBox->registerParent(_Parent);
  };
  inline void deregisterParent(Box &_Parent) {
    if (apBox) apBox->deregisterParent(_Parent);
  };

  inline bool operator < (OrderedAxisElm const &_E) const {
#ifndef NDEBUG
    assert(apBox != NULL && _E.apBox != NULL);
#endif
    double const &v1 = getVal();
    double const &v2 = _E.getVal();
    return v1!=v2 ? v1<v2 : getId()!=_E.getId() ? getId()<_E.getId() : getType()<_E.getType();
  };

  inline bool operator > (OrderedAxisElm const &_E) const {
    return ! operator < (_E);
  };

  inline bool operator == (OrderedAxisElm const &_E) const {
    return getKey() == _E.getKey();
  };


  Box::T_Key aKey;
  bool aYAxisSearch : 1;
  bool aMaxSearch : 1;
  bool aFixed : 1;
  double aVal;
  Box *apBox;


  friend ostream & operator << (ostream &, OrderedAxisElm const &);
};


ostream & operator << (ostream &_s, OrderedAxisElm const &_o) {
#ifndef NDEBUG
#if MY_DEBUG > 4
  _s
    << "{" << _o.getId() 
    << ":" << Box::typeAsStr(_o.getType())
    << ":" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? ">" : "<")
    << ":Fixed:" << _o.isFixed()
    << ":Val:" << _o.getVal()
    << ": Ref:" << _o.apBox
    << ": Ref.ParentsNum:" << _o.apBox->aParents.size()
    << "}";
#elif MY_DEBUG > 2
  _s
    << "{" << _o.getId() 
    << ":" << Box::typeAsStr(_o.getType())
    << ":" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? ">" : "<")
    << ":Fixed:" << _o.isFixed()
    << ":Val:" << _o.getVal()
    << "}";
#else
  _s
    << "{" << _o.getId() 
    << ":" << Box::typeAsStr(_o.getType())
    << ":" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? ">" : "<")
    << ":" << _o.isFixed()
    << ":(" << _o.getVal()
    << ")}";

#endif
#endif
  return _s;
}



/**
 * Comparator for std::set to order axis.
 * Matches the default one - no need to use.
 */
struct OrderedAxisElmCmp {
  inline bool operator ()(OrderedAxisElm const &_A, OrderedAxisElm const &_B) const {
    return _A < _B;
  };
};


/**
 * Hasher.
 * Needed for unordered_set.
 */
struct OrderedAxisElmHash {
  inline size_t operator ()(OrderedAxisElm const &_A) const {
    Box::T_Key key = _A.getKey();
    //return key.first << key.second;
    return hash<int>{}(key.first) ^ hash<int>{}(key.second);
  };
};


// typedef set<OrderedAxisElm, AxisCmp> T_AxisOrdered;
// typedef map<double, Box*> T_AxisOrdered;
// typedef vector<OrderedAxisElm> T_AxisOrdered;
typedef map<OrderedAxisElm, bool> T_AxisOrdered;
typedef unordered_set<OrderedAxisElm, OrderedAxisElmHash> T_AxisUnordered;
typedef T_AxisOrdered::value_type T_AxisOrderedElm;


ostream & operator << (ostream &_s, T_AxisOrderedElm const &_o) {
  return _s << _o.first;
};




/**********************************************************************************
 * ################################################################################
 * Scene region
 */


typedef map<int, BoxedVrtx *> T_BoxPartMap;
typedef pair<int, BoxedVrtx *> T_MapPair;


typedef vector<BoxedVrtx> T_BoxPartCntr;
typedef vector<BoxedEdge> T_BoxEdgeCntr;
typedef vector<BoxedHalfP> T_BoxHalfPCntr;

typedef map<int, BoxedVrtx *> T_BoxPartMap;
typedef pair<int, BoxedVrtx *> T_MapPair;




/**
 * Scene region.
 * It is assumed that all children (if any) cover the whole space of it's
 * parent, thus a parent having children does not store objects itself.
 */
class BoxedRegion : public Box
{
  public:

    typedef vector<BoxedRegion> T_Regions;
    typedef pair<Box::T_Key, bool> T_ActiveCntrElm;
    typedef unordered_map<Box::T_Key, bool, Box::KeyHash> T_ActiveCntr;


    struct BoxAdder {
      Box &aBox;
      BoxAdder(Box &_Box) : aBox(_Box) { };
      void operator ()(BoxedRegion &_R) { _R.add(aBox); };
    };

    struct BoxUpdater {
      Box &aBox;
      BoxUpdater(Box &_Box) : aBox(_Box) { };
      void operator ()(BoxedRegion &_R) { _R.update(aBox); };
    };

    struct BoxEraser {
      Box &aBox;
      BoxEraser(Box &_Box) : aBox(_Box) { };
      void operator ()(BoxedRegion &_R) { _R.erase(aBox); };
    };

    struct BoxPEraser {
      Box &aBox;
      BoxPEraser(Box &_Box) : aBox(_Box) { };
      void operator ()(BoxedRegion &_R) { _R.erasePermanently(aBox); };
    };

    struct IntersectFinder {
      bool Res = false;
      T_IntersCntr &PP, &PE, &PH;
      IntersectFinder(T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) : PP(_PP), PE(_PE), PH(_PH) { };
      void operator ()(BoxedRegion &_R) { Res |= _R.findIntersect(PP, PE, PH); };
    };

    struct RegionsCounter {
      size_t Num;
      RegionsCounter() : Num(0) {};
      void operator ()(BoxedRegion const &_R) { Num += _R.countRegions(); };
    };

    struct RegionsSplitter {
      size_t RegionsGen;
      RegionsSplitter(int _RegionsGen) : RegionsGen(_RegionsGen) {};
      void operator ()(BoxedRegion &_R) { _R.split(RegionsGen); };
    };


    static int getObjPerRegionPow2() {
      return cObjPerRegLimitPow2;
    };
    static size_t getObjPerRegion() {
      return 1L << cObjPerRegLimitPow2;
    }

    static int cObjPerRegLimitPow2;



    BoxedRegion() :
      Box(), aDevidedYAxis(true),  aAutoSplit(true), apBoxedVrtxs(NULL), apBoxedEdges(NULL), apBoxedHalfP(NULL) {
      };
    BoxedRegion(int _Id,
        T_BoxPartCntr &_BoxedVrtxs,
        T_BoxEdgeCntr &_BoxedEdges,
        T_BoxHalfPCntr &_BoxedHalfP) :
      Box(_Id, E_RegionBox), aDevidedYAxis(true), aAutoSplit(true),
      apBoxedVrtxs(&_BoxedVrtxs), apBoxedEdges(&_BoxedEdges), apBoxedHalfP(&_BoxedHalfP) {
        setInfBoundry();
        setInitialized();
        changed();
      };


    inline void clear() {
      aXAxis.clear();
      aYAxis.clear();
      aChildren.clear();
    };

    inline bool hasChildren() const {
      return aChildren.size();
    };
    inline bool isLeaf() const {
      return !hasChildren();
    };

    inline void propagateResetChange() {
      resetChange();
      for_each(aXAxis.begin(), aXAxis.end(), [](T_AxisOrderedElm &_E) { _E.second = false; });
      for_each(aXAxis.begin(), aXAxis.end(), [](T_AxisOrderedElm &_E) { _E.second = false; });
      for_each(aChildren.begin(), aChildren.end(), [](BoxedRegion &_R) { _R.propagateResetChange(); });
    };

    inline void add(Box &_B) {
      if (!intersects(_B)) {
        return;
      } else if (hasChildren()) {
        BoxAdder adder(_B);
        for_each(aChildren.begin(), aChildren.end(), adder);
      } else {
        addToOrderedAxis(_B);
        _B.registerParent(*this);
      };
    };

    inline void erase(Box &_B) {
      if (hasChildren()) {
        BoxEraser eraser(_B);
        for_each(aChildren.begin(), aChildren.end(), eraser);
      } else {
        eraseFromOrderedAxis(_B);
      };
    };
    inline void remove(Box &_B) {
      if (hasChildren()) {
        for_each(aChildren.begin(), aChildren.end(), [&](BoxedRegion &_R) { _R.remove(_B); });
      } else {
        removeOnUpdateFromOrderedAxis(_B);
      };
    };
    inline void update(Box &_B) {
      if (!intersects(_B)) {
        erasePermanently(_B);
      } else if (hasChildren()) {
        for_each(aChildren.begin(), aChildren.end(), [&](BoxedRegion &_R) { _R.update(_B); });
      } else {
        updateOrderedAxis(_B);
        _B.registerParent(*this);
      };
    };
    inline void finalizeUpdate() {
      if (hasChildren()) {
        for_each(aChildren.begin(), aChildren.end(), [](BoxedRegion &_R) { _R.finalizeUpdate();});
      } else {
        finalizeUpdateLoc();
      };
    }

    inline void erasePermanently(Box &_B) {
      if (hasChildren()) {
        BoxPEraser eraser(_B);
        for_each(aChildren.begin(), aChildren.end(), eraser);
      } else {
#ifndef NDEBUG
#if MY_DEBUG > 2
        cout << "Removing permanently:" << _B << " from:" << *this << endl;
#endif
#endif
        eraseFromOrderedAxis(_B);
        _B.deregisterParent(*this);
      };
    };

    int preSplit(int _Num, Box const &_Scene);

    inline bool split(int _RegionsGen) {
      if (!_RegionsGen)
        return false;

      _RegionsGen--;
      if (hasChildren()) {
        return splitChildren(_RegionsGen);
      } else if (splitLoc(_RegionsGen)) {
        return splitChildren(_RegionsGen);
      };

      return false;
    };

    void merge();

    inline bool findIntersect(T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) {
      if(checkAndSplit()) {
        IntersectFinder finder(_PP, _PE, _PH);
        for_each(aChildren.begin(), aChildren.end(), finder);
        return finder.Res;
      } else {
        return findIntersectLoc(_PP, _PE, _PH);
      };
    };

    inline size_t countRegions() const {
      //size_t Num=1;
      //return Num + (aChildren.size() ? aChildren.at(0).countRegions() + aChildren.at(1).countRegions() : 0);
      RegionsCounter counter;
      for_each(aChildren.begin(), aChildren.end(), counter);
      return counter.Num;
    };

    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedRegion const &);


  private:

    inline void addToOrderedAxis(Box &_O) {
      // for 'add' both axis are maintained - needed for split
      aXAxis.insert(T_AxisOrderedElm(OrderedAxisElm(false, false, _O),true));
      aXAxis.insert(T_AxisOrderedElm(OrderedAxisElm(false, true, _O),true));
      aYAxis.insert(T_AxisOrderedElm(OrderedAxisElm(true, false, _O),true));
      aYAxis.insert(T_AxisOrderedElm(OrderedAxisElm(true, true, _O),true));
    };

    inline void eraseFromOrderedAxis(Box &_O) {
      aXAxis.erase(OrderedAxisElm(false, false, _O));
      aXAxis.erase(OrderedAxisElm(false, true, _O));
      aYAxis.erase(OrderedAxisElm(true, false, _O));
      aYAxis.erase(OrderedAxisElm(true, true, _O));
    };

    inline void updateOrderedAxis(Box &_O) {
      // for update only one axis is maintained - no need for both
      aXAxis.insert(T_AxisOrderedElm(OrderedAxisElm(false, false, _O),true));
      aXAxis.insert(T_AxisOrderedElm(OrderedAxisElm(false, true, _O),true));
    };

    inline void removeOnUpdateFromOrderedAxis(Box &_O) {
      aXAxis.erase(OrderedAxisElm(false, false, _O));
      aXAxis.erase(OrderedAxisElm(false, true, _O));
    };

    inline void finalizeUpdateLoc() {
      // negative performance affect!
      //aXAxis.clear(); aYAxis.clear();
      //copy(aXAxisHash.begin(), aXAxisHash.end(), inserter(aXAxis, aXAxis.begin()));
      //copy(aYAxisHash.begin(), aYAxisHash.end(), inserter(aYAxis, aYAxis.begin()));
    };

    inline double checkAxisSpan(int _Size, T_AxisOrdered &_Axis) {
      // return (*_Axis.rbegin()).first - (*_Axis.begin()).first;
      return (*_Axis.rbegin()).first.getVal() - (*_Axis.begin()).first.getVal();
    };

    inline bool checkXSplit(int _XSize, int _YSize) {
      return checkAxisSpan(_XSize, aXAxis) > checkAxisSpan(_YSize, aYAxis);
    };

    inline T_AxisOrdered const & getMinAxis() const {
      // return aDevidedYAxis ? aYAxis : aXAxis;
      return aXAxis.size() > aYAxis.size() ? aYAxis : aXAxis;
    };
    inline T_AxisOrdered & getMinAxis() {
      // return aDevidedYAxis ? aYAxis : aXAxis;
      return aXAxis.size() > aYAxis.size() ? aYAxis : aXAxis;
    };

    inline bool shouldSplit() const {
      // min Axis contains ~ 2 * objects (min,max)
      return aXAxis.size() > 10 and aYAxis.size() > 10 && (getMinAxis().size()>>1) > getObjPerRegion();
    };

    inline int getEstObjs() const {
      return getMinAxis().size()>>1;
    };

    inline bool checkAndSplit() {
      if (isLeaf()) {
        if (aAutoSplit && shouldSplit()) {
#ifndef NDEBUG
#if MY_DEBUG > 0
          cout << __FUNCTION__
            << " auto split of " << (*this)
            << ", EstObjs=" << getEstObjs()
            << ", ObjPerReg=" << getObjPerRegion()
            << endl;
#endif
#endif
          return split(log2(getEstObjs() / getObjPerRegion()));
        };
        return false;
      } else {
        return true;
      };
    };

    inline bool splitChildren(int _RegionsGen) {
      RegionsSplitter splitter(_RegionsGen);
      for_each(aChildren.begin(), aChildren.end(), splitter);
      return true;
    };

    bool splitAxis(int, double &_Devider,
        Box &_A, Box &_B,
        T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
        T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis, 
        T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis);

    bool findAxisIntersectLoc(T_AxisOrdered &, T_IntersCntr &, T_IntersCntr &, T_IntersCntr &);
    bool findIntersectLoc(T_IntersCntr &, T_IntersCntr &, T_IntersCntr &);
    bool splitLoc(int );


  private:

    bool aDevidedYAxis;
    bool aAutoSplit;

    T_AxisOrdered aXAxis;
    T_AxisOrdered aYAxis;
    //T_AxisUnordered aXAxisHash;
    //T_AxisUnordered aYAxisHash;

    T_BoxPartCntr *apBoxedVrtxs;
    T_BoxEdgeCntr *apBoxedEdges;
    T_BoxHalfPCntr *apBoxedHalfP;

    T_Regions aChildren;

    T_IntersCntr aPP;
    T_IntersCntr aPE;
    T_IntersCntr aPH;
};



int BoxedRegion::cObjPerRegLimitPow2 = 30;


ostream & BoxedRegion::toStr(ostream &_s) const
{
  Box::toStr(_s)
    << ", Leaf:" << !hasChildren()
    << ", DevidedYAxis:" << aDevidedYAxis
    << ", ChildrenNum:" << aChildren.size();

  dumpContainer<>(g_MaxCoutNum, _s, NULL, "\n ****  XAxis:",
      aXAxis.size(), aXAxis.begin(), aXAxis.end());

  dumpContainer<>(g_MaxCoutNum, _s, NULL, "\n ****  YAxis:",
      aYAxis.size(), aYAxis.begin(), aYAxis.end());

  if (hasChildren()) {
    _s << "\n **** (" << getId() << ")->" << aChildren[0]
      << "\n **** (" << getId() << ")->" << aChildren[1];
  };

  return _s;
}


ostream & operator << (ostream &_s, BoxedRegion const &_o)
{
  return _o.toStr(_s);
}


ostream & operator << (ostream &_s, BoxedRegion::T_ActiveCntrElm const &_o)
{
  return _s << _o.first << "(Fixed:" << _o.second << ")";
}



/**
 *
 * Only for the first (inf) region before any split or registration.
 */
int BoxedRegion::preSplit(int _Num, Box const &_Scene)
{
  int Counter = 0;

  if (!isInitialized() || hasChildren()) {
    return Counter;
  };

  clear();
  int IdBase = getId();

  for(int i=0; i < 8; ++i) {
    T_Boundry Boundry = getBoundry();
    switch(i) {
      case 0:
        Boundry[1] = _Scene.xmin();
        Boundry[3] = _Scene.ymin();
        break;
      case 1:
        Boundry[1] = _Scene.xmin();
        Boundry[2] = _Scene.ymin();
        Boundry[3] = _Scene.ymax();
        break;
      case 2:
        Boundry[1] = _Scene.xmin();
        Boundry[2] = _Scene.ymax();
        break;
      case 3:
        Boundry[0] = _Scene.xmin();
        Boundry[1] = _Scene.xmax();
        Boundry[2] = _Scene.ymax();
        break;
      case 4:
        Boundry[0] = _Scene.xmax();
        Boundry[2] = _Scene.ymax();
        break;
      case 5:
        Boundry[0] = _Scene.xmax();
        Boundry[2] = _Scene.ymin();
        Boundry[3] = _Scene.ymax();
        break;
      case 6:
        Boundry[0] = _Scene.xmax();
        Boundry[3] = _Scene.ymin();
        break;
      case 7:
        Boundry[0] = _Scene.xmin();
        Boundry[1] = _Scene.xmax();
        Boundry[3] = _Scene.ymin();
        break;
    };

    IdBase++;
    aChildren.push_back(BoxedRegion(IdBase, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP));
    BoxedRegion &R = aChildren.back();
    R.registerParent(*this);
    R.setBoundry(Boundry);
    ++Counter;
  };

  int Max = _Num-1;
  Box RegionBox = _Scene;
  Box RegionBox2 = _Scene;
  RegionBox2.shrinkBy(_Num);

  for (int i=0; i < Max ; ++i) {
    RegionBox.setBoundry(RegionBox2);
    for (int j=0; j < Max ; ++j) {
      IdBase++;
      aChildren.push_back(BoxedRegion(IdBase, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP));
      BoxedRegion &R = aChildren.back();
      R.registerParent(*this);
      R.setBoundry(RegionBox);
      ++Counter;

      if (!RegionBox.safeShiftX()) {
        clear();
        return 0;
      };
    };
    if (!RegionBox2.safeShiftY()) {
      clear();
      return 0;
    };
  };

  RegionBox2.setYMax(_Scene.ymax());

  for (int i=0; i < Max; ++i) {
    IdBase++;
    aChildren.push_back(BoxedRegion(IdBase, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP));
    BoxedRegion &R = aChildren.back();
    R.registerParent(*this);
    R.setBoundry(RegionBox2);
    ++Counter;

    if (!RegionBox2.safeShiftX()) {
      clear();
      return 0;
    };
  };

  RegionBox2.setXMax(_Scene.xmax());
  IdBase++;
  aChildren.push_back(BoxedRegion(IdBase, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP));
  BoxedRegion &R = aChildren.back();
  R.registerParent(*this);
  R.setBoundry(RegionBox2);
  ++Counter;

  RegionBox2 = _Scene;
  RegionBox2.setXMin(RegionBox.xmin());
  RegionBox2.shrinkYBy(_Num);

  for (int i=0; i < Max; ++i) {
    IdBase++;
    aChildren.push_back(BoxedRegion(IdBase, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP));
    BoxedRegion &R = aChildren.back();
    R.registerParent(*this);
    R.setBoundry(RegionBox2);
    ++Counter;

    if (!RegionBox2.safeShiftY()) {
      clear();
      return 0;
    };
  };

#ifndef NDEBUG
#if MY_DEBUG > 2
  cout << __FUNCTION__
    << " SplitFactor=" << _Num
    << " SubRegions=" << Counter
    << " Scene=" << _Scene << endl;
  dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, "\n **** Children:",
      aChildren.size(), aChildren.begin(), aChildren.end(), "\n") << endl;
#endif
#endif

  return Counter;
};




/**
 * Splits region in half, either by X or Y
 *
 * @returns     True when split,
 *              False when reached limit of objects per region or no more generations requested
 */
bool BoxedRegion::splitLoc(int _RegionsGen)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif

  int Num = 2;

#ifndef NDEBUG
#if MY_DEBUG > 2
  cout << __FUNCTION__ << " RegionsGenerations:" << _RegionsGen << " for:" << *this << endl;
#endif
#endif

  if (!isInitialized())
  {
    return false; 
  }
  else if (shouldSplit())
  {
    int IdBase = getId() << 1;
    for (int i=0; i < Num; ++i) {
      BoxedRegion R(IdBase+i+1, *apBoxedVrtxs, *apBoxedEdges, *apBoxedHalfP);
      aChildren.push_back(R);
      aChildren.back().registerParent(*this);
    };

    int XSize = aXAxis.size();
    int YSize = aYAxis.size();
    int BoundryIdx = 0;
    double Devider = 0;

    BoxedRegion &A = aChildren[0];
    BoxedRegion &B = aChildren[1];
    T_Boundry ABoundry = getBoundry();
    T_Boundry BBoundry = getBoundry();

    // 
    // Split if any objects sorted
    if (XSize && YSize)
    {
      //
      // should check of the validity of the split (objects might be big and
      // cover both forked regions) indeally a new divider might be checked,
      // now just rollback
      //

      if (checkXSplit(XSize, YSize)) {
        if(!splitAxis(XSize, Devider,
              A, B,
              aXAxis, aYAxis,
              A.aXAxis, A.aYAxis,
              B.aXAxis, B.aYAxis)) {
          return false;
        };
        aDevidedYAxis = false;
      } else {
        if(!splitAxis(YSize, Devider,
              A, B,
              aYAxis, aXAxis, 
              A.aYAxis, A.aXAxis, 
              B.aYAxis, B.aXAxis)) {
          return false;
        };
        aDevidedYAxis = true;
      };

      aXAxis.clear();
      aYAxis.clear();
    };

    BoundryIdx = aDevidedYAxis << 1;

    A.aDevidedYAxis = !aDevidedYAxis;
    B.aDevidedYAxis = !aDevidedYAxis;

    ABoundry[BoundryIdx+1] = Devider;
    BBoundry[BoundryIdx] = Devider;

    A.setBoundry(ABoundry);
    B.setBoundry(BBoundry);

#ifndef NDEBUG
#ifdef MY_TIMING
    D_TIME_SEC2(time, "split Region", getKey());
#endif
#endif

    return true;
  };

  return false;
}



/**
 * Splits axis for spawned children regions.
 * Regions borders "touch", i.e. the split line belongs to both regions.
 */
bool BoxedRegion::splitAxis(int _ASize, double &_Devider,
    Box &_A, Box &_B,
    T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
    T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis,
    T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis)
{
  bool Ret = false;
  int SplitMark = _ASize >> 1;
  vector<OrderedAxisElm> SplitAxis;
  vector<OrderedAxisElm> ASplitAxis;
  vector<OrderedAxisElm> BSplitAxis;
  unordered_set<OrderedAxisElm, OrderedAxisElmHash> Active;

  SplitAxis.clear();
  SplitAxis.reserve(_SplitSrcAxis.size());
  ASplitAxis.clear();
  ASplitAxis.reserve(_SplitSrcAxis.size());
  BSplitAxis.clear();
  BSplitAxis.reserve(_SplitSrcAxis.size());

  transform(_SplitSrcAxis.begin(), _SplitSrcAxis.end(), back_inserter(SplitAxis), 
      [](T_AxisOrderedElm const &_E){ return _E.first;});
  SplitAxis.shrink_to_fit();

  // (SplitMark > 0) - checked earlier
  OrderedAxisElm const &I = SplitAxis[SplitMark-1];
  OrderedAxisElm const &J = SplitAxis[SplitMark];
  _Devider = (J.getVal() + I.getVal())/2.0;
  //_Devider = ((*SplitAxis.begin()).getVal() + (*SplitAxis.rbegin()).getVal())/2.0;

  vector<OrderedAxisElm>::iterator Iter = SplitAxis.begin();
  copy_n(Iter, SplitMark, back_inserter(ASplitAxis));
  Iter = SplitAxis.begin();
  advance(Iter, SplitMark);
  copy(Iter, SplitAxis.end(), back_inserter(BSplitAxis));

  Active.clear();
  Iter = SplitAxis.begin();
  advance(Iter, SplitMark);
  for_each(SplitAxis.begin(), Iter, [&](OrderedAxisElm &_E) { 
      if (_E.isMaxSearch())
      Active.erase(_E);
      else
      Active.insert(_E);
      });
  // add all openeded
  move(Active.begin(), Active.end(), back_inserter(BSplitAxis));

  Active.clear();
  vector<OrderedAxisElm>::reverse_iterator RIterFrom = SplitAxis.rbegin();
  vector<OrderedAxisElm>::reverse_iterator RIterTo = SplitAxis.rbegin();
  advance(RIterTo, SplitAxis.size() - SplitMark);
  for(;RIterFrom != RIterTo; ++RIterFrom) {
    if ((*RIterFrom).isMaxSearch())
      Active.insert((*RIterFrom));
    else
      Active.erase((*RIterFrom));
  };
  // add all openeded
  move(Active.begin(), Active.end(), back_inserter(ASplitAxis));

  ASplitAxis.shrink_to_fit();
  BSplitAxis.shrink_to_fit();

  if (ASplitAxis.size() < SplitAxis.size() && BSplitAxis.size() < SplitAxis.size())
  {
    sort(ASplitAxis.begin(), ASplitAxis.end());
    transform(ASplitAxis.begin(), ASplitAxis.end(), inserter(_ASplitAxis, _ASplitAxis.begin()),
        [](OrderedAxisElm const &_E) { return T_AxisOrderedElm(_E, false);});
    sort(BSplitAxis.begin(), BSplitAxis.end());
    transform(BSplitAxis.begin(), BSplitAxis.end(), inserter(_BSplitAxis, _BSplitAxis.begin()),
        [](OrderedAxisElm const &_E) { return T_AxisOrderedElm(_E, false);});

    unordered_set<Box::T_Key, Box::KeyHash> AKeys, BKeys;
    transform(ASplitAxis.begin(), ASplitAxis.end(),
        inserter(AKeys, AKeys.begin()),
        [](OrderedAxisElm const &_E) { return _E.getKey(); });
    transform(BSplitAxis.begin(), BSplitAxis.end(),
        inserter(BKeys, BKeys.begin()),
        [](OrderedAxisElm const &_E) { return _E.getKey(); });

    copy_if(_KeepSrcAxis.begin(), _KeepSrcAxis.end(),
        inserter(_AKeepAxis, _AKeepAxis.begin()),
        [&](T_AxisOrderedElm const &_E) { return AKeys.find(_E.first.getKey())!=AKeys.end(); });
    copy_if(_KeepSrcAxis.begin(), _KeepSrcAxis.end(),
        inserter(_BKeepAxis, _BKeepAxis.begin()),
        [&](T_AxisOrderedElm const &_E) { return BKeys.find(_E.first.getKey())!=BKeys.end(); });

    //for_each(ASplitAxis.begin(), ASplitAxis.end(), [&](OrderedAxisElm &_E) { if(_E.isMaxSearch()) _E.registerParent(_A); });
    //for_each(BSplitAxis.begin(), BSplitAxis.end(), [&](OrderedAxisElm &_E) { if(_E.isMaxSearch()) _E.registerParent(_B); });
    //for_each(SplitAxis.begin(), SplitAxis.end(), [&](OrderedAxisElm &_E) {_E.deregisterParent(*this); });

#ifndef NDEBUG
#if MY_DEBUG > 2
    cout << __FUNCTION__
      << " SplitSize=" << _ASize
      << " SplitMark=" << SplitMark
      << " Devider=" << _Devider << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig   SplitAxis:",
        _SplitSrcAxis.size(), _SplitSrcAxis.begin(), _SplitSrcAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Tmp    SplitAxis:",
        SplitAxis.size(), SplitAxis.begin(), SplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child ASplitAxis:",
        _ASplitAxis.size(), _ASplitAxis.begin(), _ASplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child BSplitAxis:",
        _BSplitAxis.size(), _BSplitAxis.begin(), _BSplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig   KeepAxis:",
        _KeepSrcAxis.size(), _KeepSrcAxis.begin(), _KeepSrcAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child AKeepAxis:",
        _AKeepAxis.size(), _AKeepAxis.begin(), _AKeepAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child BKeepAxis:",
        _BKeepAxis.size(), _BKeepAxis.begin(), _BKeepAxis.end()) << endl;
#endif
#endif
    Ret = true;
  }
  else
  {
    aChildren.clear();

#ifndef NDEBUG
#if MY_DEBUG > 1
    cout << __FUNCTION__
      << " Ineffective split for:" << (*this)
      << endl;
#endif
#endif
    Ret = false;
  };

  return Ret;
}



void BoxedRegion::merge()
{
}



bool BoxedRegion::findAxisIntersectLoc(
    T_AxisOrdered &_Axis,
    T_IntersCntr &_PP,
    T_IntersCntr &_PE,
    T_IntersCntr &_PH)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif

  bool bFound = false;

  struct IntersectsInserter {
    IntersectsInserter(
        vector<OrderedAxisElm>::iterator _I,
        T_IntersCntr &_PP,
        T_IntersCntr &_PE,
        T_IntersCntr &_PH,
        T_BoxPartCntr const *_pVrtxCntr,
        T_BoxEdgeCntr const *_pEdgeCntr,
        T_BoxHalfPCntr const *_pHalfPCntr) : 
      apElm(_I),
      aPP(_PP),
      aPE(_PE),
      aPH(_PH),
      apVrtxCntr(_pVrtxCntr),
      apEdgeCntr(_pEdgeCntr),
      apHalfPCntr(_pHalfPCntr)
    {
    };

    void operator ()(OrderedAxisElm const &_Elm)
    {
      if (_Elm.isMaxSearch()) return;
      if (!(apElm->getBox()->intersectsY(*_Elm.getBox()))) return;
      if ((*apElm).isFixed() && _Elm.isFixed()) return;

      switch((apElm->getType() | _Elm.getType())) {
        case E_VrtxBox: 
          {
            aPP.push_back( (apElm->getId() < _Elm.getId() ?
                  pair<int,int>(apElm->getId(), _Elm.getId()) :
                  pair<int,int>(_Elm.getId(), apElm->getId())));
          }; break;
        case E_EdgeBox:
          break;
        case 3:
          {
            int Id1, Id2;
            if (apElm->getType() == E_EdgeBox) {
              Id1 = _Elm.getId(); Id2 = apElm->getId();
            } else {
              Id1 = apElm->getId(); Id2 = _Elm.getId();
            }
            if(!apEdgeCntr->at(Id2).isMember(Id1)) {
              aPE.push_back(pair<int,int>(Id1, Id2));
            };
          }; break;
        case E_HalfPBox:
          break;
        case 5:
          {
            aPH.push_back((apElm->getType() == E_HalfPBox ?
                  pair<int,int>(_Elm.getId(), apElm->getId()) :
                  pair<int,int>(apElm->getId(), _Elm.getId())));
          }; break;
      };

#ifndef NDEBUG
#if MY_DEBUG > 4
      cout << __FUNCTION__
        << ", Akey:" << *apElm
        << ", Bkey:" << _Elm
        << ", Intersecs:" << intersects
        << ", IntersectsSize:" << (intersects ? intersects->size() : -1)
        << endl;
#endif
#endif
    };

    vector<OrderedAxisElm>::iterator apElm;
    T_IntersCntr &aPP, &aPE, &aPH;
    T_BoxPartCntr const *apVrtxCntr;
    T_BoxEdgeCntr const *apEdgeCntr;
    T_BoxHalfPCntr const *apHalfPCntr;
  };


  vector<OrderedAxisElm> WeededAxis;
  // WeededAxis.reserve(_Axis.size());
  /*copy_if(_Axis.begin(), _Axis.end(), back_inserter(WeededAxis), 
      [&](T_AxisOrdered::value_type const &_E) { 
        AxisMark += _E.isMaxSearch() ? -1 : 1;
        return AxisMark > 0 && (_E.getKey() != (*++I).getKey());
      });*/
  transform(_Axis.begin(), _Axis.end(), back_inserter(WeededAxis),
      [](T_AxisOrderedElm const &_E) { return _E.first; });
  /*
  for_each(_Axis.begin(), _Axis.end(),
      [&](T_AxisOrderedElm const &_E) { if(_E.second) WeededAxis.push_back(_E.first); });
      */


  // 
  // Looping over axis to find intersections
  // T_AxisOrdered::iterator I;
  vector<OrderedAxisElm>::iterator I;
  OrderedAxisElm AxisKey;
  IntersectsInserter Inserter(WeededAxis.end(), _PP, _PE, _PH, apBoxedVrtxs, apBoxedEdges, apBoxedHalfP);
  for (I=WeededAxis.begin(); I != WeededAxis.end(); ++I) {
    if (!(*I).isMaxSearch()) {
      if (Inserter.apElm != WeededAxis.end()) {
        for_each(I, find(I, WeededAxis.end(), AxisKey.initAxisPair(*(Inserter.apElm))), Inserter);
      };
      Inserter.apElm = I;
    } else if (Inserter.apElm != WeededAxis.end() && (*Inserter.apElm).getKey() == (*I).getKey()) {
      Inserter.apElm = WeededAxis.end();
    };
  };

  /*
  T_AxisOrdered::iterator AI = _Axis.begin();
  T_AxisOrdered::iterator AI2 = _Axis.end();
  vector< T_AxisOrdered::iterator > tmp;
  for (;AI != _Axis.end(); AI++) {
    AI2 = AI; AI++;
    if (!(AI2->first.isMaxSearch())) {
      if (AI->first.isMaxSearch() && AI2->first.getKey() == AI->first.getKey())
        continue;
      else
        tmp.push_back(AI2);
    } else {
    };
  };
  */

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC4(time, "findIntersect Axis", getKey(), _Axis.size(), WeededAxis.size());
#endif
#endif

  return bFound;
}


bool BoxedRegion::findIntersectLoc(
    T_IntersCntr &_PP,
    T_IntersCntr &_PE,
    T_IntersCntr &_PH)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time_total = clock();
#endif
#endif

  size_t PPSize = _PP.size();
  size_t PESize = _PE.size();
  size_t PHSize = _PH.size();

  if (hasChanged() && isInitialized())
  {
    aPP.clear(); aPE.clear(); aPH.clear();
    findAxisIntersectLoc(aXAxis, aPP, aPE, aPH);

#ifndef NDEBUG
#if MY_TIMING > 0
    D_TIME_SEC4(time_total, "findIntersect XY", getKey(), aXAxis.size(), aYAxis.size());
#endif
#endif


#ifndef NDEBUG
#if MY_DEBUG > 0
    cout << __FUNCTION__ 
      << " PPx:" << aPP.size()
      << " PEx:" << aPE.size()
      << " PHx:" << aPH.size()
      << " Region: " << getKey()
      << endl;
#if MY_DEBUG > 3
    cout << __FUNCTION__ 
      << "\n** Intersections for Region:" << (*this);
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** aPP:",
        aPP.size(), aPP.begin(), aPP.end());
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** aPE:",
        aPE.size(), aPE.begin(), aPE.end());
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** aPH:",
        aPH.size(), aPH.begin(), aPH.end()) << endl;
#endif
#endif
#endif

  }
  else
  {

#ifndef NDEBUG
#if MY_DEBUG > 0
    cout << __FUNCTION__
      << " No changes: " << (*this)
      << endl;
#endif
#endif

  };

  // restoring last frame collisions if no changes
  _PP.reserve(aPP.size());
  copy(aPP.begin(), aPP.end(), back_inserter(_PP));
  _PE.reserve(aPE.size());
  copy(aPE.begin(), aPE.end(), back_inserter(_PE));
  _PH.reserve(aPH.size());
  copy(aPH.begin(), aPH.end(), back_inserter(_PH));

#ifndef NDEBUG
#if MY_TIMING > 0
  D_TIME_SEC(time_total, "Copy containers");
#endif
#endif

#ifndef NDEBUG
#if MY_DEBUG > 2
  cout << __FUNCTION__ 
    << " Region: " << getKey();
  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** _PP:",
      _PP.size(), _PP.begin(), _PP.end());
  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** _PE:",
      _PE.size(), _PE.begin(), _PE.end()) << endl;
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
class BoxedScene : Box
{

  public:

    BoxedScene() : Box(), aRootRegion(0, aBoxVrtxs, aBoxEdges, aBoxHalfP) {};

    inline bool isInit() const { return aInitialized; };

    inline void init(TwoDScene const &_Scene) {
      boxObjs(_Scene);
    };
    void update(TwoDScene const &);

    inline bool findIntersect(T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) {
      return aRootRegion.findIntersect(_PP, _PE, _PH);
    };

    inline void addToRegion(Box &_Box) {
      aRootRegion.add(_Box);
    };
    inline void eraseFromRegion(Box &_Box) {
      aRootRegion.erase(_Box);
    };
    inline void updateToRegion(Box &_Box) {
      aRootRegion.update(_Box);
    };
    inline void removeFromRegion(Box &_Box) {
      aRootRegion.remove(_Box);
    };

    inline int getObjNum() const {
      return aVrtxsNum + aEdgesNum + aHalfPsNum;
    };
    inline int getEstRegionsNum() const {
      return getObjNum() >> BoxedRegion::getObjPerRegionPow2(); // ~256 per region
    };
    inline int getEstRegionsGen() const {
      int Gens = log2(getEstRegionsNum());
      return Gens > 0 ? Gens : 0;
    };
    inline size_t getRegionsNum() const {
      return aRootRegion.countRegions();
    };


    //friend ostream & operator << (ostream &_, BoxedVrtx const &_o);


  private:

    void boxVrtxs(TwoDScene const &);
    void boxEdges(TwoDScene const &);
    void boxHalfP(TwoDScene const &);
    void boxObjs(TwoDScene const &);
    
    inline void precalcBoundry(Box &_B) {
      for_each(aBoxVrtxs.begin(), aBoxVrtxs.end(), [&](BoxedVrtx const &_V) { _B.expandNoInf(_V);});
    };
    inline int preSplitRegion(int _ObjNum, Box &_SceneBoundry) {
      int Num = 2; // _ObjNum >> 4 / aRootRegion.getObjPerRegion();
      return aRootRegion.preSplit(Num, _SceneBoundry);
    };

    void updateVrtx(TwoDScene const &, VectorXs const &, BoxedVrtx &);
    void updateEdge(TwoDScene const &, BoxedEdge &);


  private:

    bool aInitialized;
    double aAvgSize;
    int aVrtxsNum;
    int aEdgesNum;
    int aHalfPsNum;

    BoxedRegion aRootRegion;

    T_BoxPartCntr aBoxVrtxs; 
    T_BoxEdgeCntr aBoxEdges; 
    T_BoxHalfPCntr aBoxHalfP;

    T_BoxPartMap aParticleMap;
};



/**
 * Boxes scene objects.
 */
void BoxedScene::boxVrtxs(TwoDScene const &_Scene)
{
  VectorXs const &X = _Scene.getX();

  aBoxVrtxs.clear();
  aBoxVrtxs.reserve(aVrtxsNum);
  for (int i=0; i < aVrtxsNum; ++i) {
    BoxedVrtx V(i, _Scene.isFixed(i), X.segment<2>(i<<1), _Scene.getRadius(i));
    aBoxVrtxs.push_back(V);
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxEdges(TwoDScene const &_Scene)
{
  aBoxEdges.clear();
  aBoxEdges.reserve(aEdgesNum);
  for (int i=0; i < aEdgesNum; ++i) {
    double R = _Scene.getEdgeRadii()[i];
    T_IntersIds Vrtxs = _Scene.getEdge(i);
    BoxedEdge E(i, &aBoxVrtxs.at(Vrtxs.first), &aBoxVrtxs.at(Vrtxs.second), R);
    E.setEdgeFixed();
    aBoxEdges.push_back(E);
    aBoxEdges.back().registerSelf();
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxHalfP(TwoDScene const &_Scene)
{
  aBoxHalfP.clear();
  aBoxHalfP.reserve(aHalfPsNum);
  for (int i=0; i < aHalfPsNum; ++i) {
    pair<VectorXs,VectorXs> const &HP = _Scene.getHalfplane(i);
    BoxedHalfP H(i, HP.first.segment<2>(0), HP.second.segment<2>(0));
    aBoxHalfP.push_back(H);
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxObjs(TwoDScene const &_Scene)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif

  VectorXs const &X = _Scene.getX();

  aAvgSize = 0.0;

  aVrtxsNum = _Scene.getNumParticles();
  aEdgesNum = _Scene.getNumEdges();
  aHalfPsNum = _Scene.getNumHalfplanes();

  Box InitSceneBoundry(0, E_RegionBox);
  boxVrtxs(_Scene);
  precalcBoundry(InitSceneBoundry);
  int PreSplit = 0; // makes things worse preSplitRegion(aVrtxsNum, InitSceneBoundry);

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC2(time, "Prespliting scenery", PreSplit);
#endif

  cout << __FUNCTION__
    << " Presplit result:" << PreSplit
    << ", RootRegion:" << aRootRegion
    << endl;
#endif

  for_each(aBoxVrtxs.begin(), aBoxVrtxs.end(), [&](BoxedVrtx &_V) { addToRegion(_V);});

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC2(time, "Boxing vertixes", aVrtxsNum);
#endif
#endif

  // Edges have to go after Vertexes!
  boxEdges(_Scene);
  for_each(aBoxEdges.begin(), aBoxEdges.end(), [&](BoxedEdge &_E) { addToRegion(_E);});

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC2(time, "Boxing edges", aEdgesNum);
#endif
#endif

  boxHalfP(_Scene);
  for_each(aBoxHalfP.begin(), aBoxHalfP.end(), [&](BoxedHalfP &_H) { addToRegion(_H);});

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC2(time, "Boxing halfplanes", aHalfPsNum);
#endif
#endif

  aRootRegion.split(getEstRegionsGen());

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC2(time, "Spliting regions with pNum", aVrtxsNum);
#endif
#endif
  aInitialized = true;


#ifndef NDEBUG

  cout << __FUNCTION__
    << ", Particles=" << aVrtxsNum
    << ", Edges=" << aEdgesNum
    << ", HalfPs=" << aHalfPsNum
    << ", AvgSize=" << aAvgSize
    << ", TotalSceneObj=" << getObjNum()
    << ", EstObjPerReg=" << getEstRegionsNum()
    << ", EstRegGen=" << getEstRegionsGen()
    << ", RegionsNum=" << getRegionsNum();

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n ****  VrtxBox:",
      aBoxVrtxs.size(), aBoxVrtxs.begin(), aBoxVrtxs.end());

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n ****  EdgeBox:",
      aBoxEdges.size(), aBoxEdges.begin(), aBoxEdges.end());

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** HalfPBox:",
      aBoxHalfP.size(), aBoxHalfP.begin(), aBoxHalfP.end());

  cout << "\n" << aRootRegion
    << endl;

#endif
}





void BoxedScene::updateVrtx(TwoDScene const &_Scene, VectorXs const &_X, BoxedVrtx &_O) {
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

    removeFromRegion(_O);
    _O.update(X, _Scene.getRadius(_O.getId()));
    updateToRegion(_O);
    _O.changed();
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

    removeFromRegion(_O);
    _O.update(_Scene.getEdgeRadii()[_O.getId()]);
    updateToRegion(_O);
    _O.changed();
  };
};



/**
 * Updates data from scene changes.
 */
void BoxedScene::update(TwoDScene const &_Scene)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif

  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  for_each(aBoxVrtxs.begin(), aBoxVrtxs.end(), [](BoxedVrtx &_o) { _o.resetChange(); });
  for_each(aBoxEdges.begin(), aBoxEdges.end(), [](BoxedEdge &_o) { _o.resetChange(); });
  aRootRegion.propagateResetChange();

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "reset change");
#endif
#endif


  T_BoxPartCntr::iterator pIter = aBoxVrtxs.begin();
  for (; pIter != aBoxVrtxs.end(); ++pIter) {
    if (!(*pIter).isFixed()) {
      updateVrtx(_Scene, X, *pIter);
    };
  };


#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "update vrtxs");
#endif
#endif


  // Edges have to go after Vertexes!
  T_BoxEdgeCntr::iterator pIterE = aBoxEdges.begin();
  for (; pIterE != aBoxEdges.end(); ++pIterE) {
    if (!(*pIter).isFixed()) {
      updateEdge(_Scene, *pIterE);
    };
  };

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "update edges");
#endif
#endif


#ifndef NDEBUG
#ifdef MY_DEBUG
#if MY_DEBUG > 1

  cout << aRootRegion << endl;

#endif
#endif
#endif
}



static BoxedScene g_Scene;




/**********************************************************************************
 * ################################################################################
 * Local module functions.
 */




static void estimateEuler(
    bool _First,
    TwoDScene const &_Scene,
    PPList &_PP,
    PEList &_PE,
    PHList &_PH)
{
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif


  if (!_First) {
    g_Scene.update(_Scene);
  };

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "update");
#endif
#endif

  T_IntersCntr PP, PE, PH;
  g_Scene.findIntersect(PP, PE, PH);

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "findIntersect");
#endif
#endif

  if (PP.size()) {
    sort(PP.begin(), PP.end());
    copy(PP.begin(), PP.end(), inserter(_PP, _PP.begin()));
  };

  if (PE.size()) {
    sort(PE.begin(), PE.end());
    copy(PE.begin(), PE.end(), inserter(_PE, _PE.begin()));
  };

  if (PH.size()) {
    sort(PH.begin(), PH.end());
    copy(PH.begin(), PH.end(), inserter(_PH, _PH.begin()));
  };

#ifndef NDEBUG
#ifdef MY_TIMING
  D_TIME_SEC(time, "copy with inserter");
#endif
#endif

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
  int HalfPsNum = _Scene.getNumHalfplanes();

  combination(ParticlesNum, _PP);
  combination2(ParticlesNum, EdgesNum, _PE);
  combination2(ParticlesNum, HalfPsNum, _PH);
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
#ifndef NDEBUG
#ifdef MY_TIMING
  clock_t time = clock();
#endif
#endif

  _First = !g_Scene.isInit();
  if (_First) {
    g_Scene.init(_Scene);

#ifndef NDEBUG
#ifdef MY_TIMING
    D_TIME_SEC(time, "Init scene");
#endif
#endif

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


static ostream & operator << (ostream &_s, T_IntersCntr const &_o) {
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

#ifndef NDEBUG
#if MY_DEBUG > 3

  if (pppairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** got PP:",
        pppairs.size(), pppairs.begin(), pppairs.end()) << endl;
  if (pepairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** got PE:",
        pepairs.size(), pepairs.begin(), pepairs.end()) << endl;
  if (phpairs.size())
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** got PH:",
        phpairs.size(), phpairs.begin(), phpairs.end()) << endl;

#endif
#endif

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
