//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 29/03/2014

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <math.h>

using namespace std;

int C[10000001], P[1000001], rq, nf, np, nd;
long long int n, q, F[50], D[10000];


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

void Divisores(long long int q)
{    int nda;long long int k;
     F[0] = 1;  nd = 1;  D[1] = 1;  nda = 1;
     for (int i=1; i<=nf; i++)
     {   if (F[i] == F[i-1]) k = k*F[i];
         else
         {    nda = nd;  k = F[i];
         }
         for (int j=1; j<=nda; j++)
             D[++nd] = k*D[j];
     }
}


int main()
{
    long long int n, cont, U, L, aux, m, m2, h;
    cin >> n;
    GeraCrivo(31323);
    GeraPrimos(31323);
    while(n--)
    {
        cin >> U >> L;

        m = 0;
        m2 = 0;
        aux = U;
        for(h = U; h <= L; h++)
        {
            Fatora(h);
            Divisores(h);
            if(nd > m)
            {
                        m = nd;
                        m2 = h;

            }

        }
        cout << "Between " << U << " and " << L << ", "  << m2 << " has a maximum of " << m << " divisors.\n";



    }


    return 0;
}
