//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 30/03/2014

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
{   unsigned int cont, ri, i;
    rq = 46341;  n = rq;  n = n*n;
    GeraCrivo(rq);
    GeraPrimos(rq);
    while(cin >> q)
    {
        nf = 0;
        if(TestaPrimo(q)) cout << (q-1)/2 << endl;
        else
        {
            Fatora(q);
            F[0] = 0;
            r = 0;
            ri = 1;
            cont = 1;
            i=nf;
            while(F[i])
            {
                if(F[i] == F[i-1])
                    cont++;
                else
                {
                    ri *= (pow(F[i],(cont)) - pow(F[i],(cont-1)));
                    cont =1;
                }
                i--;
            }

            cout << ri/2 << endl;

        }

    }
    return 0;
}
