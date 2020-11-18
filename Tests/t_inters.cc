#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>



using namespace std;


ostream & operator << (ostream &_s, pair<int,int> const &_o) {
  return _s << "{" << _o.first << "," << _o.second << "}";
};

ostream & operator << (ostream &_s, vector<pair<int,int> > const &_o) {
  //copy(_o.begin(), _o.end(), ostream_iterator<pair<int,int> >(_s, ", "));
  for(size_t i=0; i < _o.size(); ++i) {
    _s << _o.at(i);
  };

  return _s;
}

ostream & operator << (ostream &_s, vector<int> const &_o) {
  copy(_o.begin(), _o.end(), ostream_iterator<int>(_s, ", "));
  return _s;
}


/**
 * Finds 2 vectors insersection.
 * Vectors have to be sorted
 * @param _A a smaller sized vector
 * @param _B a bigger sized vector
 * @param _C destination
 */
template<typename T>
inline bool findSortedCntrIntersect(T const &_A, T const &_B, T &_C) {
  typename T::const_iterator IA;
  typename T::const_iterator IB1, IB2;
  register int OrigSizeC = _C.size();

  cout << "v1:" << _A << "\nv2:" << _B << "\nv3:" << _C << endl;
  for (IA=_A.cbegin(), IB1=_B.cbegin(); IA!=_A.cend(); ++IA) {
    IB2 = find(IB1, _B.cend(), (*IA));
    cout << "checking: " << *IA << " and " << *IB2 << endl;
    if (IB2 != _B.cend()) {
      _C.push_back((*IA));
      cout << "##### found intersection:" << (*IA) << endl;
      IB1 = IB2;
    };
  };

  cout << "v1:" << _A << "\nv2:" << _B << "\nv3:" << _C << endl;
  return _C.size() - OrigSizeC;
};


/**
 * Finds 2 vectors insersection.
 * Vectors will be sorted
 * @param _A a vector
 * @param _B a vector
 * @param _C destination
 */
template<typename T>
inline bool findCntrIntersect(T &_A, T &_B, T &_C) {
  register int SizeA, SizeB;
  SizeA = _A.size(); SizeB = _B.size();
  if(SizeA && SizeB) {
    sort(_A.begin(), _A.end());
    sort(_B.begin(), _B.end());

    if (SizeA < SizeB) {
      _C.reserve(SizeA);
      return findSortedCntrIntersect(_A, _B, _C);
    } else {
      _C.reserve(SizeB);
      return findSortedCntrIntersect(_B, _A, _C);
    };
  };
  _C.shrink_to_fit();

  return false;
}


void t1() {
  typedef pair<int,int> T_VE;
  typedef vector<T_VE> T_V;

  T_V v1 = {{2,3},{1,2},{4,9}};
  T_V v2 = {{0,1},{2,3},{1,2},{1,1}};
  T_V v3;

  findCntrIntersect(v1,v2,v3);

  cout << "V3:" << v3;
  sort(v2.begin(), v2.end());
  cout << "\nsorted V2:" << v2 << endl;
}

void t2() {
  typedef int T_VE;
  typedef vector<T_VE> T_V;

  T_V v1 = {2,5,6,8,4,3,5,7};
  T_V v2 = {9,3,2};
  T_V v3;

  findCntrIntersect<vector<int> >(v1,v2,v3);

  cout << "V3:";
  copy(v3.begin(), v3.end(), ostream_iterator<T_VE>(cout, ", "));
  cout << endl;
}


int main()
{

  t1();
  t2();


  return 0;
}
