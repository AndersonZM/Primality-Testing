//AC em 30/04/2014
//Anderson Zudio
//Victor Cracel
//Julio Cesar

//       -------------------ATENCAO!!!!!!----------------------
//       ESSE PROBLEMA ACEITA PREDIOS DE TAMANHO NEGATIVO!!!!!!
//       -------------------ATENCAO!!!!!!----------------------

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <math.h>

using namespace std;
int T[10001], S[10001], W[10001], n, m, l, kscm;

int main()
{
    int cont = 0, t, aux, decre, cre;
    cin >> t;
    while(t--)
    {
        cin >> n;
        cre = decre = 0;
        for(int i = 1; i <= n; i++ )
            cin >> S[i];

        for(int i = 1; i <= n; i++)
        {
            cin >> W[i];
            if(W[i] < 0)
            {
                for(int j = i; j <= n; j++)
                    S[j] = S[j+1];
                i--;
                n--;
            }
        }

        fill(T, T+10001, 0);
        for (int i = 1; i <= n; i++)
        {
             T[i] = W[i];
             for (int j = 1; j < i; j++)
             {
                if (S[i] > S[j])
                {
                    T[i] = max(T[i], T[j] + W[i]);
                }
             }
         }

        for(int i = 1; i <= n; i++)
        {
            cre = max(cre, T[i]);
        }


        for(int i = 1, j = n; i <= floor(n/2); i++, j--)
        {
            aux = S[j];
            S[j] = S[i];
            S[i] = aux;
            aux = W[j];
            W[j] = W[i];
            W[i] = aux;
        }

        fill(T, T+10001, 0);
        for (int i = 1; i <= n; i++)
        {
             T[i] = W[i];
             for (int j = 1; j < i; j++)
             {
                if (S[i] > S[j])
                {
                    T[i] = max(T[i], T[j] + W[i]);
                }
             }
         }

        for(int i = 1; i <= n; i++)
        {
            decre = max(decre, T[i]);
        }

        if(cre >= decre)
            cout << "Case " << ++cont << ". Increasing (" << cre << "). " << "Decreasing (" << decre << ")." << endl;
        else
            cout << "Case " << ++cont << ". Decreasing (" << decre << "). " << "Increasing (" << cre << ")." << endl;
    }
    return 0;
}
