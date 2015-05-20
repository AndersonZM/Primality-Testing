//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014


#include <iostream>
#include <math.h>
using namespace std;

int K, seq, i, cont, cont3, aux[1000100];
bool achei;

int main()
{
   while(cin >> K)
   {
            seq = 3;
            cont = 3;
            aux[3] = 1;
            achei = true;
            while(achei)
            {
                seq *= 10;
                seq += 3;
                seq = seq % K;
                if(aux[seq])
                {
                    cout << cont << " " << cont - 2 << " " << cont - 1 << endl;
                    achei = false;
                    break;
                }
                cont++;
            }


   }


return 0;
}
