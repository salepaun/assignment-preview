#include <iostream>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>
#include <ctime>
#include <vector>
#include <iterator>

#include <string.h>


using namespace std;


ostream & operator << (ostream &_s, vector<int> const &_v) {
  copy(_v.begin(), _v.end(), ostream_iterator<int>(_s, ","));
  return _s;
};



int main(int _c, char **_v) {

  double a1[] = {1.1, 1.2};
  vector<int> v1 = {1,2,3,4,5,6,7,8,9};
  int i = v1.size() >> 1;

  vector<int> v2, v3, v4, v5;

  v2.reserve(v1.size());
  v3.reserve(v1.size());
  v4.reserve(v1.size());
  v5.reserve(v1.size());

  vector<int>::iterator FI;
  vector<int>::reverse_iterator RI1, RI2;

  copy(v1.begin(), v1.end(), back_inserter(v5, v5.begin()));

  copy_n(v1.begin(), i, inserter(v2, v2.begin()));
  FI = v1.begin();
  advance(FI, i);
  copy(FI, v1.end(), inserter(v3, v3.begin()));

  RI1 = v1.rbegin();
  RI2 = v1.rbegin();
  advance(RI2, v1.size()-i);
  for(;RI1 != RI2; ++RI1) {
    v4.push_back(*RI1);
  };

  cout << "i=" << i
    << "\nv1=" << v1
    << "\nv2=" << v2
    << "\nv3=" << v3
    << "\nv4=" << v4
    << "\nv5=" << v5
    << "\n double[2] size=" << sizeof(a1) << ", sizeof(double)=" << sizeof(double)
    << endl;


  return 0;
}
