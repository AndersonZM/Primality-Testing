//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 29/03/2014

#include <iostream>
#include <math.h>
using namespace std;

int C[10000001], P[1000001], rq, nf, np, nd;
long long int n, q, F[50], D[10000], r;
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

void Fatora (long long int q)
{    nf = 0;
     for (int i=1; i<=np; i++)
     {   while((q % P[i]) == 0)
         {   F[++nf] = P[i];  q = q/P[i];
         }
         if (q == 1) break;
     }
     if (q != 1) F[++nf] = q;
}

int main()
{   int cont, ri;
    rq = 10000000;  n = rq;  n = n*n;
    GeraCrivo(rq);
    GeraPrimos(rq);
    while(1)
    {   cin >> q;
        if(!q) break;
        nf = 0;
        if(q < 0) q = -q;
        Fatora(q);
        int i;
        cont = 0;
        for (i=nf; i> 0; i--)
        {
            teste = TestaPrimo(F[i]);
            if (teste)
            {
                if(!cont)
                {
                    r = F[i];
                    ri = i;
                    cont++;
                }
                else if(F[i] != r)
                    cont++;
            }

        }
        if(cont >= 2) cout << F[ri] << endl;
        else cout << "-1" << endl;
    }
    return 0;
}







