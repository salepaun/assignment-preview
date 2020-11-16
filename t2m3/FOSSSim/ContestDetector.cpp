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
#include <cmath>

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
#define MY_DEBUG 1


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

    inline double const & xmin() const { return aBorder[0]; };
    inline double const & xmax() const { return aBorder[1]; };
    inline double const & ymin() const { return aBorder[2]; };
    inline double const & ymax() const { return aBorder[3]; };

    inline bool intersects(Box const &_B) const {
      return intersects(_B.getBoundry());
    };
    inline bool intersects(T_Boundry const &_B) const {
      bool R = intersectsAxis(_B[0], _B[1], aBorder[0], aBorder[1])
        && intersectsAxis(_B[2], _B[3], aBorder[2], aBorder[3]);
      return R;
    };
    inline bool intersectsAxis(double const &_An, double const &_Ax,
        double const &_Bn, double const &_Bx) const {
      return !(_An > _Bx || _Ax < _Bn);
    };
    inline bool intersectsAxisOld(double const &_An, double const &_Ax,
        double const &_Bn, double const &_Bx) const {
      return (_An <= _Bx && _An >= _Bn) || (_Bn <= _Ax && _Bn >= _An);
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
class BoxedObj : public Box
{

  public:

    BoxedObj() : Box() {};
    BoxedObj(int _Id) : Box(_Id, E_VrtxBox) {};
    BoxedObj(int _Id, bool _Fixed, Vector2s const &_X, double const &_R) : 
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

    friend ostream & operator << (ostream &, BoxedObj const &);


  private:

    double aX;
    double aY;
};



ostream & BoxedObj::toStr(ostream &_s) const {
  return Box::toStr(_s)
    << " (" << aX << "," << aY
    << ")";
}


ostream & operator << (ostream &_s, BoxedObj const &_o) {
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

    typedef pair<BoxedObj *, BoxedObj *> T_VrtxPair;

    BoxedEdge() : Box(), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id) : Box(_Id, E_EdgeBox, false), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id, BoxedObj *_A, BoxedObj *_B, double const &_R) : 
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

    inline BoxedObj const * a() const {
      return aVrtx.first;
    };
    inline BoxedObj const * b() const {
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
 * Boxed Objects Container supporting data.
 */



/**
 * Element for ordered axis.
 */
struct OrderedAxisElm
{
  OrderedAxisElm() : aKey(Box::T_Key(-1, E_UnknownBox)), aYAxisSearch(false), aMaxSearch(false), aFixed(false), apBox(NULL)
  { };

  OrderedAxisElm(bool _YAxis, bool _Max, Box &_B) :
    aYAxisSearch(_YAxis),
    aMaxSearch(_Max),
    apBox(&_B)
  { 
    aKey = _B.getKey();
    aFixed = _B.isFixed();
    aVal = getRawVal(*this, _B);
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
  return _s
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
  return _s
    << "{" << _o.getId() 
    << ":" << Box::typeAsStr(_o.getType())
    << ":" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? ">" : "<")
    << ":Fixed:" << _o.isFixed()
    << ":Val:" << _o.getVal()
    << "}";
#else
  return _s
    << "{" << _o.getId() 
    << ":" << Box::typeAsStr(_o.getType())
    << ":" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? ">" : "<")
    << ":" << _o.isFixed()
    << ":(" << _o.getVal()
    << ")}";

#endif
#endif
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
    return key.first << key.second;
  };
};


// typedef set<OrderedAxisElm, AxisCmp> T_AxisOrdered;
typedef set<OrderedAxisElm> T_AxisOrdered;







/**********************************************************************************
 * ################################################################################
 * Scene region
 */

class BoxedHalfP : public Box {
};

typedef map<int, BoxedObj *> T_BoxPartMap;
typedef pair<int, BoxedObj *> T_MapPair;


typedef vector<BoxedObj> T_BoxPartCntr;
typedef vector<BoxedEdge> T_BoxEdgeCntr;
typedef vector<BoxedHalfP> T_BoxHalfPCntr;

typedef map<int, BoxedObj *> T_BoxPartMap;
typedef pair<int, BoxedObj *> T_MapPair;




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
    typedef map<Box::T_Key, bool> T_ActiveCntr;


    struct BoxAdder {
      Box &aBox;
      BoxAdder(Box &_Box) : aBox(_Box) { };
      void operator ()(BoxedRegion &_R) { _R.add(aBox); };
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

    struct AxisIntersectFinder {
      bool Res = false;
      T_AxisOrdered &A;
      T_IntersCntr &PP, &PE, &PH;
      AxisIntersectFinder(T_AxisOrdered &_A, T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) : A(_A), PP(_PP), PE(_PE), PH(_PH) { };
      void operator ()(BoxedRegion &_R) { Res |= _R.findAxisIntersect(A, PP, PE, PH); };
    };

    struct RegionsCounter {
      size_t Num;
      RegionsCounter() : Num(0) {};
      void operator ()(BoxedRegion const &_R) { Num += _R.countRegions(); };
    };



    static int getObjPerRegionPow2() {
      return cObjPerRegLimitPow2;
    };
    static size_t getObjPerRegion() {
      return 1L << cObjPerRegLimitPow2;
    }

    static int cObjPerRegLimitPow2;



    BoxedRegion() :
      Box(), aDevidedYAxis(true),  apBoxedVrtxs(NULL), apBoxedEdges(NULL), apBoxedHalfP(NULL) {
      };
    BoxedRegion(int _Id,
        T_BoxPartCntr const &_BoxedVrtxs,
        T_BoxEdgeCntr const &_BoxedEdges,
        T_BoxHalfPCntr const &_BoxedHalfP) :
      Box(_Id, E_RegionBox), aDevidedYAxis(true),
      apBoxedVrtxs(&_BoxedVrtxs), apBoxedEdges(&_BoxedEdges), apBoxedHalfP(&_BoxedHalfP) {
        setInfBoundry();
        setInitialized();
        changed();
      };


    inline bool hasChildren() const {
      return aChildren.size();
    };

    inline void propagateResetChange() {
      resetChange();
      for_each(aChildren.begin(), aChildren.end(), [](BoxedRegion &_R) { _R.propagateResetChange(); });
    };

    inline void add(Box &_B) {
      if (!intersects(_B)) {
        erasePermanently(_B);
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

    bool split(int );
    void merge();

    inline bool findAxisIntersect(T_AxisOrdered &_A, T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) {
      if (hasChildren()) {
        AxisIntersectFinder finder(_A, _PP, _PE, _PH);
        for_each(aChildren.begin(), aChildren.end(), finder);
        return finder.Res;
      } else {
        return findAxisIntersectLoc(_A, _PP, _PE, _PH);
      };
    };
    inline bool findIntersect(T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH) {
      if (hasChildren()) {
        IntersectFinder finder(_PP, _PE, _PH);
        for_each(aChildren.begin(), aChildren.end(), finder);
        return finder.Res;
      } else {
        return findIntersectLoc(_PP, _PE, _PH);
      };
    };

    inline T_AxisOrdered const & getMinAxis() const {
      // return aDevidedYAxis ? aYAxis : aXAxis;
      return aXAxis.size() > aYAxis.size() ? aYAxis : aXAxis;
    };

    inline size_t countRegions() const {
      // RegionsCounter counter;
      // for_each(aChildren.begin(), aChildren.end(), counter);
      size_t Num=1;
      return Num + (aChildren.size() ? aChildren.at(0).countRegions() + aChildren.at(1).countRegions() : 0);
    };

    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedRegion const &);


  private:

    inline void addToOrderedAxis(Box &_O) {
      aXAxis.insert(OrderedAxisElm(false, false, _O));
      aXAxis.insert(OrderedAxisElm(false, true, _O));
      aYAxis.insert(OrderedAxisElm(true, false, _O));
      aYAxis.insert(OrderedAxisElm(true, true, _O));
    };

    inline void eraseFromOrderedAxis(Box &_O) {
      aXAxis.erase(OrderedAxisElm(false, false, _O));
      aXAxis.erase(OrderedAxisElm(false, true, _O));
      aYAxis.erase(OrderedAxisElm(true, false, _O));
      aYAxis.erase(OrderedAxisElm(true, true, _O));
    };

    inline double checkAxisSpan(int _Size, T_AxisOrdered &_Axis) {
      return (*_Axis.rbegin()).getVal() - (*_Axis.begin()).getVal();
    };

    inline bool checkXSplit(int _XSize, int _YSize) {
      return checkAxisSpan(_XSize, aXAxis) > checkAxisSpan(_YSize, aYAxis);
    };

    inline bool shouldSplit() const {
      // min Axis contains ~ 2 * objects (min,max)
      return aXAxis.size() > 4 and aYAxis.size() > 4 && getMinAxis().size() > getObjPerRegion();
    };

    void splitAxis(int, double &_Devider,
        Box &_A, Box &_B,
        T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
        T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis, 
        T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis);

    bool findAxisIntersectLoc(T_AxisOrdered &, T_IntersCntr &, T_IntersCntr &, T_IntersCntr &);
    bool findIntersectLoc(T_IntersCntr &, T_IntersCntr &, T_IntersCntr &);


  private:

    bool aDevidedYAxis;

    T_AxisOrdered aXAxis;
    T_AxisOrdered aYAxis;

    T_BoxPartCntr const *apBoxedVrtxs;
    T_BoxEdgeCntr const *apBoxedEdges;
    T_BoxHalfPCntr const *apBoxedHalfP;

    T_Regions aChildren;

    T_IntersCntr aPP;
    T_IntersCntr aPE;
    T_IntersCntr aPH;
};



int BoxedRegion::cObjPerRegLimitPow2 = 1;


ostream & BoxedRegion::toStr(ostream &_s) const
{
  Box::toStr(_s)
    << ", Leaf:" << !hasChildren()
    << ", DevidedYAxis:" << aDevidedYAxis;

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
 * Splits region in half, either by X or Y
 *
 * @returns     True when split,
 *              False when reached limit of objects per region or no more generations requested
 */
bool BoxedRegion::split(int _RegionsGen)
{
  int Num = 2;

#ifndef NDEBUG
#if MY_DEBUG > 2
  cout << __FUNCTION__ << " RegionsGenerations:" << _RegionsGen << " for:" << *this << endl;
#endif
#endif

  if (!_RegionsGen)
    return false;

  _RegionsGen--;
  if (hasChildren())
  {
    for_each(aChildren.begin(), aChildren.end(), [&](BoxedRegion &_R) { _R.split(_RegionsGen); });

    return true; // To improve to collect responses - for now not used feature
  }
  else if (!isInitialized())
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

    BoxedRegion &A = aChildren.at(0);
    BoxedRegion &B = aChildren.at(1);
    T_Boundry ABoundry = getBoundry();
    T_Boundry BBoundry = getBoundry();

    // 
    // Split if any objects sorted
    if (XSize && YSize)
    {
      if (checkXSplit(XSize, YSize)) {
        splitAxis(XSize, Devider,
            A, B,
            aXAxis, aYAxis,
            A.aXAxis, A.aYAxis,
            B.aXAxis, B.aYAxis);
        aDevidedYAxis = false;
      } else {
        splitAxis(YSize, Devider,
            A, B,
            aYAxis, aXAxis, 
            A.aYAxis, A.aXAxis, 
            B.aYAxis, B.aXAxis);
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

    A.split(_RegionsGen);
    B.split(_RegionsGen);

    return true;
  }
  else
  {
#ifndef NDEBUG
#if MY_DEBUG > 1

    cout << __FUNCTION__
      << " Reached split limit for " << *this 
      << ", with AxisNum=" << getMinAxis().size()
      << endl;

#endif
#endif

    return false;
  };
}



/**
 * Splits axis for spawned children regions.
 * Regions borders "touch", i.e. the split line belongs to both regions.
 */
void BoxedRegion::splitAxis(int _ASize, double &_Devider,
    Box &_A, Box &_B,
    T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
    T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis,
    T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis)
{
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
  /* in case of changing to unordered_set 
  _ASplitAxis.reserve(_SplitSrcAxis.size());
  _BSplitAxis.reserve(_SplitSrcAxis.size());
  _AKeepAxis.reserve(_KeepSrcAxis.size());
  _BKeepAxis.reserve(_KeepSrcAxis.size());*/
  copy(_KeepSrcAxis.begin(), _KeepSrcAxis.end(), inserter(_AKeepAxis, _AKeepAxis.begin()));
  move(_KeepSrcAxis.begin(), _KeepSrcAxis.end(), inserter(_BKeepAxis, _BKeepAxis.begin()));
  move(_SplitSrcAxis.begin(), _SplitSrcAxis.end(), inserter(SplitAxis, SplitAxis.begin()));
  SplitAxis.shrink_to_fit();

  // (SplitMark > 0) - checked earlier
  OrderedAxisElm &I = SplitAxis[SplitMark-1];
  OrderedAxisElm &J = SplitAxis[SplitMark];
  _Devider = abs(J.getVal() - I.getVal())/2.0;

  vector<OrderedAxisElm>::iterator Iter = SplitAxis.begin();
  copy_n(Iter, SplitMark, inserter(ASplitAxis, ASplitAxis.begin()));
  Iter = SplitAxis.begin();
  advance(Iter, SplitMark);
  copy(Iter, SplitAxis.end(), inserter(BSplitAxis, BSplitAxis.begin()));

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
  move(Active.begin(), Active.end(), inserter(BSplitAxis, BSplitAxis.begin()));

  Active.clear();
  vector<OrderedAxisElm>::reverse_iterator RIterFrom = SplitAxis.rbegin();
  vector<OrderedAxisElm>::reverse_iterator RIterTo = SplitAxis.rbegin();
  advance(RIterTo, SplitAxis.size() - SplitMark);
  for(;RIterFrom != RIterTo; ++RIterFrom) {
    if ((*RIterFrom).isMaxSearch())
      Active.erase((*RIterFrom));
    else
      Active.insert((*RIterFrom));
  };
  // add all openeded
  move(Active.begin(), Active.end(), inserter(ASplitAxis, ASplitAxis.begin()));

  ASplitAxis.shrink_to_fit();
  BSplitAxis.shrink_to_fit();
  //for_each(ASplitAxis.begin(), ASplitAxis.end(), [&](OrderedAxisElm &_E) { if(_E.isMaxSearch()) _E.registerParent(_A); });
  //for_each(BSplitAxis.begin(), BSplitAxis.end(), [&](OrderedAxisElm &_E) { if(_E.isMaxSearch()) _E.registerParent(_B); });
  //for_each(SplitAxis.begin(), SplitAxis.end(), [&](OrderedAxisElm &_E) {_E.deregisterParent(*this); });

  sort(ASplitAxis.begin(), ASplitAxis.end());
  move(ASplitAxis.begin(), ASplitAxis.end(), inserter(_ASplitAxis, _ASplitAxis.begin()));
  sort(BSplitAxis.begin(), BSplitAxis.end());
  move(BSplitAxis.begin(), BSplitAxis.end(), inserter(_BSplitAxis, _BSplitAxis.begin()));

#ifndef NDEBUG
#if MY_DEBUG > 1
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig  SplitAxis:",
        _SplitSrcAxis.size(), _SplitSrcAxis.begin(), _SplitSrcAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Tmp   SplitAxis:",
        SplitAxis.size(), SplitAxis.begin(), SplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child SplitAxis:",
        _ASplitAxis.size(), _ASplitAxis.begin(), _ASplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child SplitAxis:",
        _BSplitAxis.size(), _BSplitAxis.begin(), _BSplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig  KeepAxis:",
        _KeepSrcAxis.size(), _KeepSrcAxis.begin(), _KeepSrcAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child KeepAxis:",
        _AKeepAxis.size(), _AKeepAxis.begin(), _AKeepAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child KeepAxis:",
        _BKeepAxis.size(), _BKeepAxis.begin(), _BKeepAxis.end()) << endl;
#endif
#endif
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
  bool bFound = false;
  T_ActiveCntr Active;

  struct IntersectsInserter {
    IntersectsInserter(T_IntersCntr &_PP, T_IntersCntr &_PE, T_IntersCntr &_PH,
        T_BoxEdgeCntr const *_pEdgeCntr) : 
      aPP(_PP),
      aPE(_PE),
      aPH(_PH),
      apEdgeCntr(_pEdgeCntr) {
      };

    void operator ()(T_ActiveCntrElm _Elm)
    {
      T_IntersCntr *intersects = NULL;
      int Id1, Id2;
      if (!aElm.second && !_Elm.second) {
        switch((aElm.first.second | _Elm.first.second)) {
          case E_VrtxBox: 
            {
              if (aElm.first.first < _Elm.first.first) {
                Id1 = aElm.first.first; Id2 = _Elm.first.first;
              } else {
                Id1 = _Elm.first.first; Id2 = aElm.first.first;
              };
              intersects = &aPP;
            }; break;
          case E_EdgeBox:
          case 3:
            {
              if (aElm.first.second == E_EdgeBox) {
                Id1 = _Elm.first.first; Id2 = aElm.first.first;
              } else {
                Id1 = aElm.first.first; Id2 = _Elm.first.first;
              }
              intersects = apEdgeCntr ? apEdgeCntr->at(Id2).isMember(Id1) ?  NULL : &aPE : &aPE;
            }; break;
          case E_HalfPBox:
          case 5:
            {
              if (aElm.first.second == E_HalfPBox) {
                Id1 = _Elm.first.first; Id2 = aElm.first.first;
              } else {
                Id1 = aElm.first.first; Id2 = _Elm.first.first;
              }
              // intersects = aEdgeCntr[Id2].isMember(Id1) ?  NULL : &aPE;
              intersects = &aPH; break;
            };
        };

#ifndef NDEBUG
#if MY_DEBUG > 2
        cout << __FUNCTION__
          << ", Akey:" << aElm
          << ", Bkey:" << _Elm
          << ", Atype|Btype:" << (aElm.first.second|_Elm.first.second)
          << ", Intersecs:" << intersects
          << endl;
#endif
#endif
        if(intersects) {
          // intersects->insert(pair<int,int>(Id1, Id2));
          intersects->push_back(pair<int,int>(Id1, Id2));
        };
      };
    };

    T_ActiveCntrElm aElm;
    T_IntersCntr &aPP, &aPE, &aPH;
    T_BoxEdgeCntr const *apEdgeCntr;
  };


#ifndef NDEBUG
#if MY_DEBUG > 2
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Axis:",
        _Axis.size(), _Axis.begin(), _Axis.end()) << endl;
#endif
#endif


  // 
  // Looping over axis to find intersections
  IntersectsInserter Inserter(_PP, _PE, _PH, apBoxedEdges);
  for (T_AxisOrdered::const_iterator I=_Axis.begin(); I != _Axis.end(); ++I) {
    Box::T_Key key = (*I).getKey();

#ifndef NDEBUG
#if MY_DEBUG > 2
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Active:",
        Active.size(), Active.begin(), Active.end()) << endl;
#endif
#endif

    if ((*I).isMaxSearch()) {
      Active.erase(key);
      continue;
    };

    if (Active.empty()) {
      Active.insert(T_ActiveCntrElm(key, (*I).isFixed()));
      continue;
    };

    T_ActiveCntrElm Elm(key, (*I).isFixed());

    Inserter.aElm = Elm;
    for_each(Active.begin(), Active.end(), Inserter);

    Active.insert(Elm);

    bFound = true;
  };

  return bFound;
}



bool BoxedRegion::findIntersectLoc(
    T_IntersCntr &_PP,
    T_IntersCntr &_PE,
    T_IntersCntr &_PH)
{
  size_t PPSize = _PP.size();
  size_t PESize = _PE.size();
  size_t PHSize = _PH.size();

  if (hasChanged()) {
    T_IntersCntr XInterPP, XInterPE, XInterPH, YInterPP, YInterPE, YInterPH;
    if (findAxisIntersect(aXAxis, XInterPP, XInterPE, XInterPH)) {
      // Improvement: only search found candidates !!!
      if (findAxisIntersect(aYAxis, YInterPP, YInterPE, YInterPH)) {
        if (XInterPP.size() && YInterPP.size()){
          sort(XInterPP.begin(), XInterPP.end()); sort(YInterPP.begin(), YInterPP.end());
          _PP.reserve(min(XInterPP.size(), YInterPP.size()));
          set_intersection(
              XInterPP.begin(), XInterPP.end(),
              YInterPP.begin(), YInterPP.end(),
              inserter(_PP, _PP.begin()));
        };

        if (XInterPE.size() && YInterPE.size()){
          sort(XInterPE.begin(), XInterPE.end()); sort(YInterPE.begin(), YInterPE.end());
          _PE.reserve(min(XInterPE.size(), YInterPE.size()));
          set_intersection(
              XInterPE.begin(), XInterPE.end(),
              YInterPE.begin(), YInterPE.end(),
              inserter(_PE, _PE.begin()));
        };

        if (XInterPH.size() && YInterPH.size()){
          sort(XInterPH.begin(), XInterPH.end()); sort(YInterPH.begin(), YInterPH.end());
          _PH.reserve(min(XInterPH.size(), YInterPH.size()));
          set_intersection(
              XInterPH.begin(), XInterPH.end(),
              YInterPH.begin(), YInterPH.end(),
              inserter(_PH, _PH.begin()));
        };
      };
    };

#ifndef NDEBUG
#if MY_DEBUG > 0
    cout << __FUNCTION__ 
      << " PPx:" << XInterPP.size()
      << " PEx:" << XInterPE.size()
      << " PHx:" << XInterPH.size()
      << " PPy:" << YInterPP.size()
      << " PEy:" << YInterPE.size()
      << " PHy:" << YInterPH.size()
      << " Changes: " << (*this)
      << endl;
#if MY_DEBUG > 3
    cout << __FUNCTION__ 
      << "** Intersections for Region:" << (*this);
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** XPP:",
        XInterPP.size(), XInterPP.begin(), XInterPP.end());
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** YPP:",
        YInterPP.size(), YInterPP.begin(), YInterPP.end());
    dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n ****  PP:",
        _PP.size(), _PP.begin(), _PP.end()) << endl;

#endif
#endif
#endif

    // saving for next frame
    aPP.reserve(_PP.size());
    copy(_PP.begin(), _PP.end(), inserter(aPP, aPP.begin()));
    aPE.reserve(_PE.size());
    copy(_PE.begin(), _PE.end(), inserter(aPE, aPE.begin()));
    aPH.reserve(_PH.size());
    copy(_PH.begin(), _PH.end(), inserter(aPH, aPH.begin()));
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

    // restoring last frame collisions if no changes
    _PP.reserve(aPP.size());
    copy(aPP.begin(), aPP.end(), inserter(_PP, _PP.begin()));
    _PE.reserve(aPE.size());
    copy(aPE.begin(), aPE.end(), inserter(_PE, _PE.begin()));
    _PH.reserve(aPH.size());
    copy(aPH.begin(), aPH.end(), inserter(_PH, _PH.begin()));
  };

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

    BoxedScene() : Box(), aRootRegion(0, aBoxVrtsx, aBoxEdges, aBoxHalfP) {};

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

    inline int getObjNum() const {
      return aVrtxsNum + aEdgesNum + aHalfPNum;
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
    int aVrtxsNum;
    int aEdgesNum;
    int aHalfPNum;

    BoxedRegion aRootRegion;

    T_BoxPartCntr aBoxVrtsx; 
    T_BoxEdgeCntr aBoxEdges; 
    T_BoxHalfPCntr aBoxHalfP;

    T_BoxPartMap aParticleMap;
};



/**
 * Boxes scene objects.
 */
void BoxedScene::boxVrtxs(TwoDScene const &_Scene)
{
  aVrtxsNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  aBoxVrtsx.clear();
  aBoxVrtsx.reserve(aVrtxsNum);
  for (int i=0; i < aVrtxsNum; ++i) {
    BoxedObj V(i, _Scene.isFixed(i), X.segment<2>(i<<1), _Scene.getRadius(i));
    addToRegion(V);
    aBoxVrtsx.push_back(V);
  };
}



/**
 * Boxes scene objects.
 */
void BoxedScene::boxEdges(TwoDScene const &_Scene)
{
  aEdgesNum = _Scene.getNumEdges();

  aBoxEdges.clear();
  aBoxEdges.reserve(aEdgesNum);
  for (int i=0; i < aEdgesNum; ++i) {
    double R = _Scene.getEdgeRadii()[i];
    T_IntersIds Vrtxs = _Scene.getEdge(i);
    BoxedEdge E(i, &aBoxVrtsx.at(Vrtxs.first), &aBoxVrtsx.at(Vrtxs.second), R);
    E.setEdgeFixed();
    addToRegion(E);
    aBoxEdges.push_back(E);
    aBoxEdges.back().registerSelf();
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

  aVrtxsNum = _Scene.getNumParticles();
  aEdgesNum = _Scene.getNumEdges();
  aHalfPNum = _Scene.getNumHalfplanes();

  boxVrtxs(_Scene);
  // Edges have to go after Vertexes!
  boxEdges(_Scene);

  aRootRegion.split(getEstRegionsGen());
  aInitialized = true;


#ifndef NDEBUG

  cout << __FUNCTION__
    << ", Particles=" << ParticlesNum
    << ", PBoxSize=" << aBoxVrtsx.size()
    << ", AvgSize=" << aAvgSize
    << ", TotalSceneObj=" << getObjNum()
    << ", EstObjPerReg=" << getEstRegionsNum()
    << ", EstRegGen=" << getEstRegionsGen()
    << ", RegionsNum=" << getRegionsNum();

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** VrtxBox:",
      aBoxVrtsx.size(), aBoxVrtsx.begin(), aBoxVrtsx.end());

  dumpContainer<>(g_MaxCoutNum, cout, NULL, "\n **** EdgeBox:",
      aBoxEdges.size(), aBoxEdges.begin(), aBoxEdges.end());

  cout << "\n" << aRootRegion
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

    eraseFromRegion(_O);
    _O.update(_Scene.getEdgeRadii()[_O.getId()]);
    addToRegion(_O);
    _O.changed();
  };
};



/**
 * Updates data from scene changes.
 */
void BoxedScene::update(TwoDScene const &_Scene)
{
  int ParticlesNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  for_each(aBoxVrtsx.begin(), aBoxVrtsx.end(), [](BoxedObj &_o) { _o.resetChange(); });
  for_each(aBoxEdges.begin(), aBoxEdges.end(), [](BoxedEdge &_o) { _o.resetChange(); });
  aRootRegion.propagateResetChange();

  T_BoxPartCntr::iterator pIter = aBoxVrtsx.begin();
  for (; pIter != aBoxVrtsx.end(); ++pIter) {
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
  T_IntersCntr PP, PE, PH;
  g_Scene.findIntersect(PP, PE, PH);

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
  _First = !g_Scene.isInit();
  if (_First) {
    g_Scene.init(_Scene);
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
