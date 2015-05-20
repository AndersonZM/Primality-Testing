//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014
#include <iostream>
#include <math.h>
using namespace std;
int a, b, g;

int C[10000001], P[1000001], F[50], D[1000], n, q, nf, np, nd;
bool teste;

void GeraCrivo (int n)
{    int i, rq, t, d;
     for (i=1; i<=n; i++)   C[i] = i;
     for (i=2; i<=n; i+=2)  C[i] = 2;
     rq = sqrt(n);
     for (i=3; i<=rq; i+=2)
     {   if (C[i] == i)
         {   t = i*i;  d = i+i;
             while(t<=n)
             {   if (C[t] == t) C[t] = i;
                 t = t+d;
             }
         }
     }
}

void GeraPrimos (int n)
{    np=0;
     for (int i=2; i<=n; i++)
         if (C[i]==i)  P[++np]= i;
}

bool TestaPrimo (int q)
{    return (C[q] == q);
}

int PotMod(int a, int b, int n){
    long long int m, c, d, rb[101];  int i, j;
    m= b;   i= 101;
    while (m > 0){
        rb[--i]= m % 2;   m= m /2;
    }
    c= 0; d= 1;
    for(j=i; j<= 100; j++){
        d= (d*d) % n;    c= 2*c;
        if (rb[j] == 1){
            d= (d*a) % n;   c++;
        }
        //cout << rb[j] << " " << c << " " << d <<endl;
    }
    return (int) d;
}

int main(){
    bool av = true;
    GeraCrivo(65000);
    GeraPrimos(65000);
    while(cin >> g, g != 0){

            if(TestaPrimo(g) == true)
            {
                cout << g << " is normal." << endl;
                continue;
            }


            av = true;
            for(int i = 2 ;i < g; i++)
            {
                if( PotMod(i, g, g) != i)
                {
                     cout << g << " is normal." << endl;
                     av = false;
                     break;
                }
            }
            if(av)
                cout << "The number " << g << " is a Carmichael number." << endl;


    }
    return 0;
}
