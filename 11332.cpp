//Anderson Zudio
//AC em 3/7/2014

#include <iostream>
using namespace std;

int main(){
  int n, x;

  while (cin >> n){
    if (n == 0)
      break;

    while(n > 9){
      x = n; n = 0;

      while (x) {
        n += x%10;
        x /= 10;
      }
    }

    cout << n << endl;
  }
  return 0;
}
