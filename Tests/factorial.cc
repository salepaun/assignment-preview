#include <iostream>


unsigned int factorial(int n) {
  int F=1;
  for (int i=2; i<=n; ++i) F*= i;
  return F;
}




int main(int _C, char **_V) {

  if (_C > 1) {
    for (int i=1; i < _C; ++i) {
      int f = atoi(*(_V+i));
      printf("%d! = %d\n", f, factorial(f));
    }
  }

  return 0;
}
