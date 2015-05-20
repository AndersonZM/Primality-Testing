//Anderson Zudio
//AC em 3/7/2014

#include <iostream>

using namespace std;

long long int x, y;

int main(){
    int t;
    cin >> t;
    while(t--){
        cin >> x >> y;
        if(x > y) cout << '>' << endl;
        else if(x < y) cout << '<' << endl;
        else cout << '=' << endl;
    }
    return 0;
}
