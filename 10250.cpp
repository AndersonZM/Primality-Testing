//Anderson Zudio de Moraes 
//Victor Cracel Messner 
//Julio Cesar Cavalcanti Carvalho
//AC em 11/04/2014 

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#define QZERO 0.0000000001
#define ABS(a,b) ((a) >= (b)) ? a-b : b-a
#define MOD(a) ((a) >= 0) ? a : -a
using namespace std;

int main(int argc, char *argv[]) {
    cout.setf(ios_base::fixed);
    cout.precision(10);
    double x1,y1,x2,y2;
    double ax1,ax2,ay1,ay2;
    double cXD, cYD, cX, cY;
    
    while (cin >> x1 >> y1 >> x2 >> y2) {
          if (((ABS(x1,x2)) <= QZERO) && ((ABS(y1,y2)) <= QZERO))
             cout << "Impossible." << endl;
          else {
               cX = (x1+x2)/2;
               cY = (y1+y2)/2;         
               cXD = (cX - x1);
               cYD = (cY - y1);
               
               ax1 = cX + cYD;
               ay1 = cY - cXD; 
               ax2 = cX - cYD;
               ay2 = cY + cXD;
              cout << ax1 << " " << ay1 << " " << ax2 << " " << ay2 << endl;
          }
    }

    return 0;
}
