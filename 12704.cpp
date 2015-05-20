//AC em 24/05/2014
//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
#include <iostream>
#include <cmath>
#include <math.h>
#define PI acos(-1)

using namespace std;




int main() {
    double t, x, y, r, dist;
    cout.setf(ios_base::fixed);
    cout.precision(2);
    cin >> t;
    while(t--) {
        cin >> x >> y >> r;
        dist = sqrt(x*x + y*y);
        cout << r - dist << " " << r + dist << endl;
    }
}
