#include <iostream>
#include <limits>
#include <cmath>


using namespace std;



int main() {
  double A = 1.1;
  double B = numeric_limits<double>::infinity();

  cout << "A=" << A << ", B=" << B << endl;
  cout << "is inf(A):" << (numeric_limits<double>::infinity() == A) \
    << ", is inf(B):" << (numeric_limits<double>::infinity() == B) \
    << ", is inf(A):" << isinf(A) \
    << ", is inf(B):" << isinf(B) \
    << endl;
}
