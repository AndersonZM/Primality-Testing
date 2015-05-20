//Anderson Zudio
//AC em 08/05/2014

#include <iostream>
#include <string.h>
#include <algorithm>

using namespace std;

int main()
{
    int t, S, n, m, K, i, j, k, l, Satual, Patual, g, h, M[305][305];
    char liberdade[305];
    while(cin >> n >> m >> K, n)
    {
        for(i = 1; i <= n; i++)
        {
            cin >> liberdade;
            for(j = 0; j < m; j++)
            {
                if(liberdade[j] == '.')
                    M[i][j+1] = 1;
                else
                    M[i][j+1] = 0;

                if(j) M[i][j+1] += M[i][j];
            }
        }

        S = n*m;
        for(j = 1; j <= m; j++)
        {
            for(k = j; k <= m; k++)
            {
                Patual = 0;
                for(i = 1, l = 1; i <= n; i++)
                {
                    Patual += M[i][k];
                    if(j > 1) Patual -= M[i][j-1];
                    if(Patual >= K)
                    {
                        while(Patual >= K && l < i)
                        {
                            Satual = (i - l + 1)*(k - j + 1);
                            S = min(S,Satual);
                            if(j > 1) Patual += M[l][j-1];
                            Patual -= M[l++][k];
                        }
                        if(Patual >= K)
                        {
                            Satual = (i - l + 1)*(k - j + 1);
                            S = min(S,Satual);
                        }
                    }
                }
            }
        }
        cout << S << endl;
    }
    return 0;
}
