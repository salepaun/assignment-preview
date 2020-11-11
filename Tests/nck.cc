#include <iostream>


using namespace std;



unsigned int factorial(int n) {
  int F=1;
  for (int i=2; i<=n; ++i) F*= i;
  return F;
}



/**
 * Calculates combinations k of n.
 */
unsigned int inline nCk(int n, int k)
{
  return factorial(n)/factorial(k)/factorial(n-k);
}




int main(int _C, char **_V) {

  if (_C > 2) {
    for (int i=1; i < _C; i+=2) {
      int n = atoi(*(_V+i));
      int k = atoi(*(_V+1+i));
      printf("%d C %d = %d\n", n, k, nCk(n,k));
    }
  }

  return 0;
}
