//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar Cavalcanti Carvalho
//AC em 30/05/2014
 
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <utility>
 
using namespace std;
 
int main() {
int T,I=0,Dx,Dy;
double a,b,c,d,t,h;
std::pair <int,int> A,B;
std::pair <double,double> C,D;
cout.setf(ios::fixed);
cout.precision(8);
    cin >> T;
    while(I < T) { cin >> A.first >> A.second >> B.first >> B.second >> b >> c >> d;
        Dx = (B.first - A.first);
        Dy = (B.second - A.second);
        a = sqrt((Dx*Dx) + (Dy*Dy));
        t = (a - c + (d * d - b * b) / (a - c)) / 2;
        h = sqrt(d * d - t * t);
 
        C.first = A.first + (Dx/a)*(t+c) - (Dy/a)*h;
        C.second = A.second + (Dy/a)*(t+c) + (Dx/a)*h;
        D.first = A.first + (Dx/a)*t - (Dy/a)*h;
        D.second = A.second + (Dy/a)*t + (Dx/a)*h;
 
        I++;
        cout.precision(8);
        cout << "Case " << I << ":" << endl << C.first << " " << C.second << " " << D.first << " " << D.second << endl;
    }
 
    return 0;
}