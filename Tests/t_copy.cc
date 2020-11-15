#include <iostream>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>
#include <ctime>

#include <string.h>


using namespace std;


#define SIZE 100000

typedef array<double, SIZE> T_ARR;


void test_copy(double (&_A)[SIZE], T_ARR &_B) {
  std::copy(_B.begin(), _B.end(), std::begin(_A));
  std::copy(std::begin(_A), std::end(_A), _B.begin());
}


void test_memcpy(double (&_A)[SIZE], T_ARR &_B) {
  memcpy(_B.data(), _A, _B.size()*sizeof(double));
  memcpy(_A, _B.data(), SIZE*sizeof(double));
}


int main(int _c, char **_v) {

  int LoopCount = 10000;
  double A[SIZE];
  array<double, SIZE> B;
  clock_t ExTimeA, ExTimeB;

  B.fill(M_PI);

  if (_c > 1) {
    LoopCount = atoi(_v[1]);
  };

  ExTimeA = clock();
  for (int i=0; i<LoopCount; ++i) {
    test_copy(A, B);
  };
  ExTimeA = clock() - ExTimeA;

  ExTimeB = clock();
  for (int i=0; i<LoopCount; ++i) {
    test_memcpy(A, B);
  };
  ExTimeB = clock() - ExTimeB;


  cout << "B size:" << B.size() << " tested copy vs memcpy on " << SIZE << " elements arrays with " << LoopCount << " loop count:"
    << "\n ** copy()   = " << (float)ExTimeA/CLOCKS_PER_SEC << " sec"
    << "\n ** memcpy() = " << (float)ExTimeB/CLOCKS_PER_SEC << " sec"
    << endl;

  return 0;
}
