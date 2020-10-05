#include "eigen3/Eigen/Dense"
#include <iostream>

using namespace std;
using namespace Eigen;

int main() {
  Matrix2f M; M << 2,0,0,2;
  VectorXf A(4); A.setZero();
  Vector2f B(1,2);
  Matrix2f N = M.inverse();

  cout << "M:\n" << M << "\ndetM=" << M.determinant() << endl;
  cout << "inversed M:\n" << N << endl;
  cout << "M*M^-1:\n" << N * M << endl;

  cout << "A:" << A << ", B:" << B << endl;
  A[0] = B[0]; A[1] = B[1];
  cout << "A:" << A << ", B:" << B << endl;

  M << 1,0,0,1;
  B << 0,10;

  cout << "M * B =" << (M * B) << endl;
  cout << "M * B * B / 2 = " << 0.5 * (M * B).dot(B) << endl;

  return 0;
}
