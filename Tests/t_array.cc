#include <array>
#include <iostream>
#include <iterator>

#include <string.h>

using namespace std;


ostream & operator << (ostream &_s, array<double, 4> const &_o) {
  copy(_o.begin(), _o.end(), ostream_iterator<double>(_s, ","));
  return _s;
}



int main() {

  double A1[4] = {1,2,3,4};
  array<double, 4> A2 = {3,4,5,6}, A3;

  memcpy(A3.data(), A1, 4*sizeof(double));

  cout << "A1:" << A1
    << ", A2:" << A2
    << ", A3:" << A3
    << endl;


  return 0;
};
