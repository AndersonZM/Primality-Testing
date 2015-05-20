//Anderson Zudio
//AC em 3/7/2014

#include <iostream>

using namespace std;

int x, j = -1, k = 2147483646;

int main(){
    int t, n;
    cin >> t;
    while(t--){
        cin >> n;
        j = -1; k = 2147483646;
        for(int i = 1; i <= n; i++){
            cin >> x;
            j = max(j, x);
            k = min(k, x);
        }
        cout << (j - k)*2 << endl;
    }
    return 0;
}
