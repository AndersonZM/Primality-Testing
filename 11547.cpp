//Anderson Zudio
//AC em 3/7/2014

#include <iostream>

using namespace std;

int main(){
    int t, x;
    cin >> t;
    while(t--){
            cin >> x;
            x *= 567; x /= 9; x += 7492; x *= 235; x /= 47; x -= 498;
            if(x < 0) x = -x;
            cout << (x/10)%10 << endl;
    }
    return 0;
}
