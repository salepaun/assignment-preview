#include <iostream>
#include <set>

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




/**
 * Creates combination '2' of n elements.
 * Clears the set.
 */
static size_t combination(int n, set<pair<int, int> > &_Set)
{
  //unsigned int PairsNum = nCk(n, 2);
  //_Set.resize(PairsNum);
  _Set.clear();
  for(int i=0; i < n; ++i)
    for(int j=i+1; j < n; ++j)
      _Set.insert(pair<int,int>(i,j));

  return _Set.size();
}


ostream & operator << (ostream &_s, pair<int,int> const &_o) {
  return _s << "(" << _o.first << "," << _o.second << ")";
};


int main(int _C, char **_V) {

  typedef set<pair<int,int> > TPairs;
  TPairs Pairs;

  if (_C > 1) {
    for (int i=1; i < _C; ++i) {
      int n = atoi(*(_V+i));
      printf("2 combination of %d = %d\n", n, combination(n, Pairs));
      for (TPairs::const_iterator I=Pairs.begin(); I != Pairs.end(); ++I) {
        cout << (*I);
      };
    };
  };

  cout << endl;

  return 0;
}
