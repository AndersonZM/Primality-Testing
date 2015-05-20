//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar Cavalcanti Carvalho
//AC em 24/05/2014
 
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#define QZERO 0.0000000001
#define ABS(a,b) ((a) >= (b)) ? a-b : b-a
#define PI acos(-1)
 
using namespace std;
   
int main() {
char city[10];
int x1, y1, z1, x2, y2 , z2;
int Px,Py,Pz, rad;
int rs, a, b , c;
double T1, T2;
cout.setf(ios::fixed);
    while(cin >> city >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> rs) {
        double delta, pctg;
        double A,B,C;
        pctg = 0.0;
        while (rs--) {
            cin >> Px >> Py >> Pz >> rad;
            a = (x2-x1);
            b = (y2-y1);
            c = (z2-z1);
            A = (a*a) + (b*b) + (c*c);
            B = 2*((a*(x1-Px)) + (b*(y1-Py)) + (c*(z1-Pz)));
            C = (x1*(x1-(2*Px))) + (y1*(y1-2*Py)) + (z1*(z1-2*Pz)) + (Px*Px) + (Py*Py) + (Pz*Pz) - (rad*rad);
           delta = (B*B) - 4*A*C;
           if (delta > 0) {
              delta = sqrt(delta);
              T1 = (-B + delta)/(2*A);
              T1 = min(max(T1,0.0),1.0);
              T2 = (-B - delta)/(2*A);
              T2 = min(max(T2,0.0),1.0);
              pctg += 100*((T1>T2) ? (T1-T2) : (T2-T1));
              //cout.precision(4);
           }
        }
        cout.precision(2);
        cout << city << endl << pctg << endl;
    }
   
    return 0;
}