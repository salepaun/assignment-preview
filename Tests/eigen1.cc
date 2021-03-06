/**
 * set makeprg=g++\ -O0\ -I/usr/include/eigen3\ -o\ %<\ %
 */



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



Vector2s cross2mv(Matrix2s const &_A, Vector2s const &_B) {
  return _B[1] * _A.col(0) - _B[0] * _A.col(1);
}


Vector2s cross2vm(Vector2s const &_A, Matrix2s const &_B) {
  return _A[0] * _B.col(1) - _A[1] * _B.col(0);
}




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

  Vector2s An = Vector2s::Random();
  Vector2s Bn = Vector2s::Random();
  Vector2s Nn = An - Bn;
  Nn.normalize();
  cout << "A:" << An.transpose() << ", B:" << Bn.transpose() \
    << ", (A - B).normalized =" << Nn << endl;

  
  Matrix2s I = Matrix2s::Identity();
  MatrixXs E(2,4);
  E.block<2,2>(0,0) = I;
  E.block<2,2>(0,2) = I;

  cout << "E:\n" << E << endl;

  cout << "Testing norm, squared norm\n"
    << " An=" << An.transpose()
    << " An.An=" << An.dot(An)
    << " An.norm=" << An.norm()
    << " An.sqNorm=" << An.squaredNorm()
    << endl;

  cout << "Testing cross product:\n"
    << " An x Bn = " << 0 // An.cross(Bn)
    << endl;

  An << 1,2;
  cout << "I.col(0)=" << I.col(0) << endl;
  cout << "I.row(0)=" << I.row(0) << endl;

  cout << "X=" << A.dot(Vector2s(I.col(0))) << endl;


  cout << "Testing Identity cross, dot:"
    << "\n I.cross(A)=" << cross2mv(I, A).transpose()
    << "\n A.cross(I)=" << cross2vm(A, I).transpose()
    << endl;

  return 0;
}
