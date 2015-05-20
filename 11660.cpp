//Anderson Zudio Moraes
//Victor Cracel Messner
//AC em 25/3/2014

#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;
int j, i;
string x[1010];

void geraSeq()
{
    int k, l = 0, w, cont;
    char atual, aux[1010];
    for(k = 1; k < i;k++)
    {
        w = 0;
        x[k] = "";
        for(l = 0; l < x[k-1].size() && l < 1001;)
        {
            atual = x[k-1][l];
            cont = 0;
            while(x[k-1][l] == atual && l < 1001)
            {
                cont++;
                l++;
            }
            sprintf(aux, "%d", cont);
            x[k] += aux;
            x[k] += atual;
        }
    }

}

int main()
{

   while(true)
   {

        cin >> x[0] >> i >> j;
        if(x[0][0] == '0' && !i && !j) break;
        geraSeq();
        cout << x[i-1][j-1] << endl;
   }

   return 0;
}
