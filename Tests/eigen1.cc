#include <iostream>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Dense>

typedef double scalar;

typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;

typedef Eigen::Matrix<scalar, 2, 2> Matrix2s;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

using namespace std;
using namespace Eigen;

int main() {
  Matrix2s M; M << 2,0,0,2;
  VectorXs A(4); A.setZero();
  Vector2s B(1,2);
  Matrix2s N = M.inverse();
  Matrix2s R = Matrix2s().Random();
  Matrix2s const & Rc = R;
  Matrix2s const & Mc = M;

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

  Vector2s Rv = M * R * B;
  cout << "10 * M * R * B = " << 10 * Rv.transpose() << endl;

  Vector2s Rvc = Mc * Rc * B;
  cout << "10 * Mc * Rc * B = " << 10 * Rvc.transpose() << endl;
  return 0;
}
