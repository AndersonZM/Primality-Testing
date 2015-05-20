//AC em 24/05/2014
//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <math.h>
#define PI acos(-1)
#define EARTH 6371.01

using namespace std;


int main() {
    int t;
    double c, f, dist, d1, d2, c1, c2, a, b, e, d;
    char X, Y, Xc, Yc;
    cin >> t;
    cout.setf(ios::fixed);
    cout.precision(2);
    while(t--){
        cin >> a >> b >> c >> X >> d >> e >> f >> Y;
        d1 = a + b/60 + c/3600;
        d1 *= PI/180;
        d2 = d + e/60 + f/3600;
        d2 *= PI/180;
        if(X == 'S') d1 = -d1;
        if(Y == 'W') d2 = -d2;
        cin >> a >> b >> c >> X >> d >> e >> f >> Y;
        c1 = a + b/60 + c/3600;
        c1 *= PI/180;
        c2 = d + e/60 + f/3600;
        c2 *= PI/180;
        if(X == 'S') c1 = -c1;
        if(Y == 'W') c2 = -c2;

        dist = cos(d1) * cos(c1) * cos(d2 - c2) + sin(d1) * sin(c1);

        if (dist < -1) dist = -1;
        else if (dist > 1) dist = 1;
        cout << acos(dist) * EARTH << endl;
    }
    return 0;
}
