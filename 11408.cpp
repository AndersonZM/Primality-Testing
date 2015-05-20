//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 30/03/2014

#include <iostream>
#include <math.h>
using namespace std;

int C[10000001], P[1000001], rq, nf, np, ndp, nd, q, k;
long long int n, F[50], D[10000] ;
bool teste, DP[5000001];

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

void Fatora (int q)
{    nf = 0;
     while(q != 1)
     {
         if(F[nf] != C[q]) F[++nf] = C[q];  q = q/C[q];
     }
}

void DePrimos(int a, int b)
{
    int soma;
    ndp = 0;
    for(int i = a; i <= b; i++)
    {
        Fatora(i);
        soma = 0;
        for(int j = 1; j <= nf; j++)
        {
            soma += F[j];
        }
        DP[++ndp] = TestaPrimo(soma);
    }
}

int main()
{
    int cont;
    n = 5000000;
    GeraCrivo(n);
    GeraPrimos(n);
    DePrimos(1,5000000);

    while(cin >> q, q != 0)
    {
        cin >> k;
        cont = 0;
        for(int i = q; i <= k; i++)
            if(DP[i]) cont++;
        cout << cont << endl;
    }
    return 0;
}







