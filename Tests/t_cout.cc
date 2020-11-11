/**
 */

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <array>

using namespace std;


typedef pair<int,int> T1;



ostream & operator << (ostream &_s, vector<int> const &_o) {
  copy(_o.begin(), _o.end(), ostream_iterator<int>(_s, " "));
  return _s;
};

ostream & operator << (ostream &_s, T1 const &_o) {
  return _s << _o.first << "," << _o.second;
};


struct Astr {
  int Id;
  bool bVar1 : 1;
  inline int id() const { return Id; };
  inline bool var1() const { return bVar1; };

  friend ostream & operator << (ostream &_s, Astr const &_o);
};

ostream & operator << (ostream &_s, Astr const &_o) {
  _s << _o.id() << ", " << _o.var1();
  return _s;
};


int main() {
  vector<pair<int,int> > cntr = { pair<int,int>(1,2) };
  vector<int> cntr2 = { 1,2,3,4,5 };
  T1 p1(1,2), p2(2,4);

  cout << cntr2 << endl;
  cout << p1 << p2 << endl;

  for(auto e: cntr) {
    cout << e << " ";
  }

  double nums[4] = {1.1, 1.2, 1.3, 1.4};
  // array<double, 4> arr(nums);

  Astr str;
  cout << str;


  /// !!! copy(cntr.begin(), cntr.end(), ostream_iterator<const T1 &>(cout, " "));
  cout << endl;
}
