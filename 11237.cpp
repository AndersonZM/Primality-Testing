//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014

#include <iostream>
#include <string.h>
#include <stdio.h>

using namespace std;

    long long int i, j, k, n, c, casas[10000100], soma, resto, osu[1000100], f, g;

int main()
{


    while(true)
    {
        cin >> c >> n;
        if((c == 0) && (n == 0)) break;

        soma = 0;
        memset(osu,0,1000100);

        for(i = 1; i <= n; i++)
        {
            cin >> casas[i];
        }

        for(i = 1; i <= n; i++)
        {
                soma += casas[i];
                resto = soma % c;
                if(!resto)
                {
                    j = 1;
                    f = i;
                    break;
                }
                else if (!(osu[resto]))
                    osu[resto]=i;
                else
                {
                    j = osu[resto] +1;
                    f = i;
                    break;
                }
        }

        for(g = j; g < f; g++)
        {
            cout << g << " ";
        }
        cout << g << endl;

    }
}
