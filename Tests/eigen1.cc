#include "eigen3/Eigen/Dense"
#include <iostream>

using namespace std;
using namespace Eigen;

int main() {
  Matrix2f M; M << 1,2,3,4;

  cout << "M:\n" << M << "\ndetM=" << M.determinant() << endl;
  cout << "inversed M:\n" << M.inverse() << endl;

  return 0;
}
