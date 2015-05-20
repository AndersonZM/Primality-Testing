//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014
#include <iostream>
#include <math.h>

using namespace std;
long long int a, p;

int C[10000001], P[1000001], rq, nf, np, nd;
long long int n, q, F[50], D[10000];
bool teste;

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
void GeraPrimos (int rq)
{    np=0;
     for (int i=2; i<=rq; i++)
         if (C[i]==i)  P[++np]= i;
}

bool TestaPrimo (long long int q)
{    bool primo = true;
     for (int i=1; i<=np; i++)
     {   if ((q % P[i]) == 0)
         {   if (q != P[i]) primo = false;
             break;
         }
     }
     return (primo);
}

int main(){

    rq = 31623;  n = rq;  n = n*n;
    GeraCrivo(rq);
    GeraPrimos(rq);

    while(true){
        cin >> p >> a;
        if (a == 0 && p == 0) break;

        if(PotMod(a, p, p) == a)
        {
            if(TestaPrimo(p))
            cout << "no" << endl;
            else cout << "yes" << endl;
        }
        else cout << "no" << endl;
    }
    return 0;
}
