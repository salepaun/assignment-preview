#include <iostream>
#include <stdio.h>

using namespace std;

int main() {
  double A1 = 10.0e-9;
  double A2 = 10.0-9;

  cout << "A1:" << A1 << endl;
  cout << "A2:" << A2 << endl;

  printf("A1: %.10lf\n", A1);
  printf("A2: %.10lf\n", A2);

  return 0;
}
