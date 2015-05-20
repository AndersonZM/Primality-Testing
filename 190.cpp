//Anderson Zudio
//Victor Messner
//AC 28/08/2014

#include <iostream>
#include <cmath>

typedef struct ponto{
    double x, y;
}ponto;

using namespace std;

int main(){
    cout.setf(ios::fixed);
    cout.precision(3);
    ponto p1, p2, p3;
    double a, b, c, k, j, m, l;
    char sinal1, sinal2, sinal3;
    while(cin >> p1.x >> p1.y >> p2.x >> p2.y >> p3.x >> p3.y){
        a = p1.x*p1.x + p1.y*p1.y;
        b = p2.x*p2.x + p2.y*p2.y;
        c = p3.x*p3.x + p3.y*p3.y;
        k= a* p2.y + c * p1.y + b * p3.y - p2.y * c - p1.y * b - p3.y * a;
        j = p2.x * a + p1.x * c + p3.x * b - p2.x * c - p1.x * b  - p3.x * a;
        l = a* p2.x * p3.y + p1.x * p2.y * c + p1.y * b * p3.x - p1.y* p2.x * c - p1.x * b * p3.y - a * p2.y * p3.x;
        m = p1.x * p2.y + p1.y * p3.x + p2.x * p3.y - p2.y * p3.x - p1.y * p2.x - p1.x * p3.y;
        double n = k/m; 
        double g = j/m;
        double u = l/m;
        if(n < 0){
            sinal1 = '+';
            n = -n;
        } else sinal1 = '-';
        if(g < 0){
            sinal2 = '-';
            g = -g;
        } else sinal2 = '+';
        double z;
        if(u > 0){
            sinal3 = '-';
            z = u;
        } else{
        	sinal3 = '+';
        	z = -u;
        } 
        
        cout << "(x " << sinal1 << " " << n/2 << ")^2 + (y " << sinal2 << " " << g/2 << ")^2 = " << sqrt((n/2)*(n/2) + (g/2)*(g/2) + u) << "^2" << endl;
        cout << "x^2 + y^2 " << sinal1 << " " << n << "x " << sinal2 << " " << g << "y " << sinal3 << " " << z << " = 0" << endl << endl;
    }
    return 0;
}