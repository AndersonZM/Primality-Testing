//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014

#include <iostream>
using namespace std;

int n, p;  long long int C[1001][1001];

void TrianguloPascal(int n)
{   int i, p;
    C[0][0] = 1; for (p=1; p<=1000; p++) C[0][p]=0;
    for (i=1; i<=n; i++)
    {   C[i][0]=1;
        for (p=1; p<=n; p++) C[i][p]=C[i-1][p]+C[i-1][p-1];
    }
}

int main()
{    TrianguloPascal(1000);
     while (1)
    {   cin >> n >> p;
        if (!(n&&p)) break;
        cout << n << " things taken " << p << " at a time is " << C[n][p] << " exactly."<<endl;
    }
}
