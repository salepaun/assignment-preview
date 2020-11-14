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


// For now: 0,1,2,3
#define MY_DEBUG 3


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


int g_MaxCoutNum = 10;



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

    static char const * typeAsStr(BoxType _Type) {
      switch (_Type) {
        case E_VrtxBox: return "V";
        case E_EdgeBox: return "E";
        case E_HalfPBox: return "H";
        case E_RegionBox: return "R";
        default: return "x";
      };
    };

    Box() : aId(-1), aType(E_UnknownBox), aChanged(false), aInitialized(false) {};
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

    inline bool intersects(Box const &_B) const {
      return intersects(_B.getBoundry());
    };
    inline bool intersects(T_Boundry const &_B) const {
      bool R = intersectsAxis(_B[0], _B[1], aBorder[0], aBorder[1])
        && intersectsAxis(_B[2], _B[3], aBorder[2], aBorder[3]);
      cout << __FUNCTION__ << " A=" << aBorder << ", B=" << _B << " R=" << R << endl;
      return R;
    };
    inline bool intersectsAxis(double const &_An, double const &_Ax,
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
      memcpy(aBorder.data(), _B, 4*sizeof(double));
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
      if (aParents.find(_Parent.getKey()) == aParents.end()) {
        aParents.insert(T_ParentElm(_Parent.getKey(), &_Parent));
        _Parent.changed();
      };
    };
    inline void deregisterParent(Box &_Parent) {
      if (aParents.find(_Parent.getKey()) != aParents.end()) {
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
          [](T_ParentElm const &_p) { _p.second->changed(); });
    };


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, Box const &);


  public:

    int aId;
    BoxType aType;
    bool aChanged;
    bool aInitialized;
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
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
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
      Box(_Id, E_VrtxBox), aFixed(_Fixed), aX(_X[0]), aY(_X[1]) {
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


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedObj const &);


  private:

    bool aFixed;
    double aX;
    double aY;
};



ostream & BoxedObj::toStr(ostream &_s) const {
  //copy(_o.aBorder.begin(), _o.aBorder.end(), ostream_iterator<double>(_s, ","));
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
    BoxedEdge(int _Id) : Box(_Id, E_EdgeBox), aVrtx(NULL, NULL) {};
    BoxedEdge(int _Id, BoxedObj *_A, BoxedObj *_B, double const &_R) : 
      Box(_Id, E_EdgeBox), aVrtx(_A, _B) {
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

    bool aFixed;
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
  OrderedAxisElm() : aYAxisSearch(false), aMaxSearch(false), aBoxObjPtr(NULL)
  { };

  OrderedAxisElm(bool _YAxis, bool _Max, Box *_ObjPtr) :
    aYAxisSearch(_YAxis),
    aMaxSearch(_Max),
    aBoxObjPtr(_ObjPtr)
  { };

  inline bool isXSearch() const { return aYAxisSearch; };
  inline bool isMaxSearch() const { return aMaxSearch; };

  inline int getId() const {
    return aBoxObjPtr ? aBoxObjPtr->getId() : -1;
  };
  inline BoxType getType() const {
    return aBoxObjPtr ? aBoxObjPtr->getType() : E_UnknownBox;
  };
  inline Box::T_Key getKey() const {
    return aBoxObjPtr ? aBoxObjPtr->getKey() : Box::T_Key(-1, E_UnknownBox);
  };
  inline double const & getVal() const {
    return getVal(*this);
  };
  inline double const & getVal(OrderedAxisElm const &_E) const {
    int i = (_E.aYAxisSearch << 1) + _E.aMaxSearch;
    return _E.aBoxObjPtr ? _E.aBoxObjPtr->aBorder[i] : aDummyVal;
  };

  inline bool operator < (OrderedAxisElm const &_E) const {
#ifndef NDEBUG
    assert(aBoxObjPtr != NULL && _E.aBoxObjPtr != NULL);
#endif
    double const &v1 = getVal();
    double const &v2 = getVal(_E);
    return v1!=v2 ? v1<v2 : getId()!=_E.getId() ? getId()<_E.getId() : getType()<_E.getType();
  };

  inline bool operator > (OrderedAxisElm const &_E) const {
    return ! operator < (_E);
  };

  inline bool operator == (OrderedAxisElm const &_E) const {
    return getKey() == _E.getKey();
  };


  double aDummyVal = 0;
  bool aYAxisSearch : 1;
  bool aMaxSearch : 1;
  Box *aBoxObjPtr;


  friend ostream & operator << (ostream &, OrderedAxisElm const &);
};


ostream & operator << (ostream &_s, OrderedAxisElm const &_o) {
  return _s
    << "{" << (_o.aYAxisSearch ? "Y" : "X")
    << ":" << (_o.aMaxSearch ? "Max" : "Min")
    << ", Id:" << _o.getId() << ":" << Box::typeAsStr(_o.getType())
    << "}";
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


typedef map<int, BoxedObj *> T_BoxPartMap;
typedef pair<int, BoxedObj *> T_MapPair;


typedef vector<BoxedObj> T_BoxPartCntr;
typedef vector<BoxedEdge> T_BoxEdgeCntr;

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
      T_IntersSet &PP, &PE, &PH;
      IntersectFinder(T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) : PP(_PP), PE(_PE), PH(_PH) { };
      void operator ()(BoxedRegion &_R) { Res |= _R.findIntersect(PP, PE, PH); };
    };

    struct AxisIntersectFinder {
      bool Res = false;
      T_AxisOrdered &A;
      T_IntersSet &PP, &PE, &PH;
      AxisIntersectFinder(T_AxisOrdered &_A, T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) : A(_A), PP(_PP), PE(_PE), PH(_PH) { };
      void operator ()(BoxedRegion &_R) { Res |= _R.findAxisIntersect(A, PP, PE, PH); };
    };

    BoxedRegion() : Box(), aDevidedYAxis(true),  apBoxedParts(NULL), apBoxedEdges(NULL) {};
    BoxedRegion(int _Id, T_BoxPartCntr const *_BoxedParts, T_BoxEdgeCntr const *_BoxedEdges) : 
      Box(_Id, E_RegionBox), aDevidedYAxis(true), apBoxedParts(_BoxedParts), apBoxedEdges(_BoxedEdges) {
        setInfBoundry();
        setInitialized();
      };


    inline bool hasChildren() const {
      return aChildren.size();
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

    void split();
    void merge();

    inline bool findAxisIntersect(T_AxisOrdered &_A, T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) {
      if (hasChildren()) {
        AxisIntersectFinder finder(_A, _PP, _PE, _PH);
        for_each(aChildren.begin(), aChildren.end(), finder);
        return finder.Res;
      } else {
        return findAxisIntersectLoc(_A, _PP, _PE, _PH);
      };
    };
    inline bool findIntersect(T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) {
      if (hasChildren()) {
        IntersectFinder finder(_PP, _PE, _PH);
        for_each(aChildren.begin(), aChildren.end(), finder);
        return finder.Res;
      } else {
        return findIntersectLoc(_PP, _PE, _PH);
      };
    };


    ostream & toStr(ostream &) const;

    friend ostream & operator << (ostream &, BoxedRegion const &);


  private:

    inline void addToOrderedAxis(Box &_O) {
      aXAxis.insert(OrderedAxisElm(false, false, &_O));
      aXAxis.insert(OrderedAxisElm(false, true, &_O));
      aYAxis.insert(OrderedAxisElm(true, false, &_O));
      aYAxis.insert(OrderedAxisElm(true, true, &_O));
      cout << *this << endl;
    };

    inline void eraseFromOrderedAxis(Box &_O) {
      aXAxis.erase(OrderedAxisElm(false, false, &_O));
      aXAxis.erase(OrderedAxisElm(false, true, &_O));
      aYAxis.erase(OrderedAxisElm(true, false, &_O));
      aYAxis.erase(OrderedAxisElm(true, true, &_O));
    };

    void splitAxis(int, double &_Devider,
        T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
        T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis, 
        T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis);

    bool findAxisIntersectLoc(T_AxisOrdered &, T_IntersSet &, T_IntersSet &, T_IntersSet &);
    bool findIntersectLoc(T_IntersSet &, T_IntersSet &, T_IntersSet &);


  private:

    bool aDevidedYAxis;

    T_AxisOrdered aXAxis;
    T_AxisOrdered aYAxis;

    T_BoxPartCntr const *apBoxedParts;
    T_BoxEdgeCntr const *apBoxedEdges;

    T_Regions aChildren;
};




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
    cout << "\n **** (" << getId() << ")->" << aChildren[0]
      << "\n **** (" << getId() << ")->" << aChildren[1];
  };

  return _s;
}


ostream & operator << (ostream &_s, BoxedRegion const &_o)
{
  return _o.toStr(_s);
}



void BoxedRegion::split()
{
  int Num = 2;

#ifndef NDEBUG
#if MY_DEBUG > 2
  cout << __FUNCTION__ << " for:" << *this << endl;
#endif
#endif

  if (hasChildren())
  {
    for_each(aChildren.begin(), aChildren.end(), [Num](BoxedRegion &_R) { _R.split(); });
  }
  else
  {
    aChildren.resize(Num);
    T_Regions::iterator I=aChildren.begin();
    int IdBase = getId() << 1;
    for (int i=0; i < Num; ++i, ++I) {
      *I = BoxedRegion(IdBase + i, apBoxedParts, apBoxedEdges);
      (*I).registerParent(*this);
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
    if (XSize || YSize)
    {
      if (XSize > YSize) {
        splitAxis(XSize, Devider,
            aXAxis, aYAxis,
            A.aXAxis, A.aYAxis,
            B.aXAxis, B.aYAxis);
        aDevidedYAxis = false;
      } else {
        splitAxis(YSize, Devider,
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
  };
}



void BoxedRegion::splitAxis(int _ASize, double &_Devider,
    T_AxisOrdered &_SplitSrcAxis, T_AxisOrdered &_KeepSrcAxis,
    T_AxisOrdered &_ASplitAxis, T_AxisOrdered &_AKeepAxis,
    T_AxisOrdered &_BSplitAxis, T_AxisOrdered &_BKeepAxis)
{
  int SplitMark = _ASize >> 1;
  vector<OrderedAxisElm> SplitAxis;

  copy(_KeepSrcAxis.begin(), _KeepSrcAxis.end(), inserter(_AKeepAxis, _AKeepAxis.begin()));
  copy(_KeepSrcAxis.begin(), _KeepSrcAxis.end(), inserter(_BKeepAxis, _BKeepAxis.begin()));
  copy(_SplitSrcAxis.begin(), _SplitSrcAxis.begin(), SplitAxis.begin());
  unordered_set<OrderedAxisElm, OrderedAxisElmHash> Active;

  //copy_n(_SplitSrcAxis.begin(), SplitMark, inserter(_ASplitAxis, _ASplitAxis.begin()));
  //copy(advance(_SplitSrcAxis.begin(), SplitMark), _SplitSrcAxis.end(), inserter(_BSplitAxis, _BSplitAxis.begin()));
  vector<OrderedAxisElm>::iterator I = SplitAxis.begin();
  for(int i; I != SplitAxis.end(); ++I, ++i) {
    if (i <= SplitMark) {
      if ((*I).isMaxSearch())
        Active.erase((*I));
      else
        Active.insert((*I));

      _ASplitAxis.insert((*I));
    };

    if (i == SplitMark) {
      // insert to second region all started not finished boxes
      copy(Active.begin(), Active.end(), inserter(_BSplitAxis, _BSplitAxis.begin()));
      _BSplitAxis.insert((*I));

      _Devider = (*I).getVal();
    };

    if (i > SplitMark) {
      if ((*I).isMaxSearch() && Active.find((*I)) != Active.end()) {
        _ASplitAxis.insert(*I);
      };
      _BSplitAxis.insert((*I));
    };
  };

#ifndef NDEBUG
#if MY_DEBUG > 1
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig  SplitAxis:",
        _SplitSrcAxis.size(), _SplitSrcAxis.begin(), _SplitSrcAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child SplitAxis:",
        _ASplitAxis.size(), _ASplitAxis.begin(), _ASplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Child SplitAxis:",
        _BSplitAxis.size(), _BSplitAxis.begin(), _BSplitAxis.end()) << endl;
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Orig  KeepAxis:",
        _SplitSrcAxis.size(), _SplitSrcAxis.begin(), _SplitSrcAxis.end()) << endl;
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
        case E_VrtxBox: 
          {
            if (aKey.first < _Key.first) {
              Id1 = aKey.first; Id2 = _Key.first;
            } else {
              Id1 = _Key.first; Id2 = aKey.first;
            };
            intersects = &aPP;
          }; break;
        case E_EdgeBox:
        case 3:
          {
            if (aKey.second == E_EdgeBox) {
              Id1 = _Key.first; Id2 = aKey.first;
            } else {
              Id1 = aKey.first; Id2 = _Key.first;
            }
            intersects = apEdgeCntr ? apEdgeCntr->at(Id2).isMember(Id1) ?  NULL : &aPE : &aPE;
          }; break;
        case E_HalfPBox:
        case 5:
          {
            if (aKey.second == E_HalfPBox) {
              Id1 = _Key.first; Id2 = aKey.first;
            } else {
              Id1 = aKey.first; Id2 = _Key.first;
            }
            // intersects = aEdgeCntr[Id2].isMember(Id1) ?  NULL : &aPE;
            intersects = &aPH; break;
          };
      };

#ifndef NDEBUG
#if MY_DEBUG > 2
      cout << __FUNCTION__
        << ", Akey:" << aKey
        << ", Bkey:" << _Key
        << ", Atype|Btype:" << (aKey.second|_Key.second)
        << ", Intersecs:" << intersects
        << endl;
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


#ifndef NDEBUG
#if MY_DEBUG > 2
    dumpContainer<>(g_MaxCoutNum, cout, __FUNCTION__, " **** Axis:",
        _Axis.size(), _Axis.begin(), _Axis.end()) << endl;
#endif
#endif


  // 
  // Looping over axis to find intersections
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



bool BoxedRegion::findIntersectLoc(
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
class BoxedScene : Box
{

  public:

    BoxedScene() : Box(), aRootRegion(0, &aBoxParts, &aBoxEdges) {};

    inline bool isInit() const { return aInitialized; };

    inline void init(TwoDScene const &_Scene) {
      boxObjs(_Scene);
    };
    void update(TwoDScene const &);

    inline bool findIntersect(T_IntersSet &_PP, T_IntersSet &_PE, T_IntersSet &_PH) {
      return aRootRegion.findIntersect(_PP, _PE, _PH);
    };

    inline void addToRegion(Box &_Box) {
      aRootRegion.add(_Box);
    };
    inline void eraseFromRegion(Box &_Box) {
      aRootRegion.erase(_Box);
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

    T_BoxPartCntr aBoxParts; 
    T_BoxEdgeCntr aBoxEdges; 

    T_BoxPartMap aParticleMap;
};



/**
 * Boxes scene objects.
 */
void BoxedScene::boxVrtxs(TwoDScene const &_Scene)
{
  aVrtxsNum = _Scene.getNumParticles();
  VectorXs const &X = _Scene.getX();

  aBoxParts.resize(aVrtxsNum);
  T_BoxPartCntr::iterator pIter = aBoxParts.begin();
  for (int i=0; i < aVrtxsNum; ++i) {
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
  aEdgesNum = _Scene.getNumEdges();

  aBoxEdges.resize(aEdgesNum);
  T_BoxEdgeCntr::iterator pIter = aBoxEdges.begin();
  for (int i=0; i < aEdgesNum; ++i) {
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

  aRootRegion.split();
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

  for_each(aBoxParts.begin(), aBoxParts.end(), [](BoxedObj &_o) { _o.resetChange(); });
  for_each(aBoxEdges.begin(), aBoxEdges.end(), [](BoxedEdge &_o) { _o.resetChange(); });

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
