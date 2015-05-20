//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar Cavalcanti Carvalho
//AC em 30/05/2014
 
#include <cstdlib>
#include <iostream>
#include <math.h>
using namespace std;
 
struct ponto {
    double x,y;
    ponto(double x,double y) : x(x), y(y) {}
    ponto() {}
};
 
int main() {
int N,I=0,K,J=0,T,i,j,p,k;
double a,b, dist, percorrido, dleft;
ponto P[101], Tree;
double D[101];
cout.setf(ios::fixed);
cout.precision(2);
    cin >> N;
    while(I < N) {
        cin >> K >> T;
        I++;
        dist = J = 0;
        while (J < K) P[J++] = (cin >> a >> b, ponto(a,b));
 
        for (i= 1; i < K ; i++) {
            a = P[i].x-P[(i-1)].x;
            b = P[i].y-P[(i-1)].y;
            D[i] = sqrt(a*a + b*b);
            dist += D[i];
        }
        dist /= (T-1);
 
        cout << "Road #" << I << ":" << endl;
        cout << P[0].x << " " << P[0].y << endl;
            p=1;
            percorrido =  0;
            dleft = dist;
            k=1;
                while (p<K) {
                    if (D[p] >= dleft) {
                        j = ceil((D[p]-dleft)/dist);
                        percorrido = dleft;
                        for (i = 0; i < j ; i++) {
                            Tree.x = P[(p-1)].x + dleft/D[p]*(P[p].x-P[(p-1)].x);
                            Tree.y = P[(p-1)].y + dleft/D[p]*(P[p].y-P[(p-1)].y);
                            cout << Tree.x << " " << Tree.y << endl;
                            k++;
                            dleft += dist;
                        }
                        dleft = dist - ((D[p]-percorrido) - ((j-1)*dist));
                    }
                    else
                        dleft = dleft - D[p];
                p++;
                }
           if (k < T)
                cout << P[(K-1)].x << " " << P[(K-1)].y << endl;
    cout << endl;
    }
    return 0;
}