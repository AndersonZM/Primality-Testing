//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014

#include <iostream>

using namespace std;

long long int V[801][2000];
int v, aux, aux2, i , j, c;
bool ver;




void GeraSubFat()
{

V[1][2000] = 0;
V[2][2000] = 1;

    for(i=3;i<=800;i++)
    {
        aux2 = 0;
        for(j=2000;j>=1;j--)
        {

            aux = V[i-1][j] + V[i-2][j] + aux2;
            V[i][j] = aux % 10;
            aux2 = aux / 10 ;
        }

        aux2 = 0;
        for(j=2000;j>=1;j--)
        {
            aux = (V[i][j])*(i-1) + aux2;
            V[i][j] = aux % 10;
            aux2 = aux/10;
        }
    }
}





int main()
{
    GeraSubFat();
    while(cin >> c)
    {
        if(c == -1) break;
        ver = 0;
        if(c == 1 ) cout << "0" << endl;
        else
        {
            for(i=1;i<= 2000;i++)
            {
                if(V[c][i] != 0) ver = 1;
                if(ver == 1) cout << V[c][i];
            }
            cout << endl;

        }
    }
return 0;
}
