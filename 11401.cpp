//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014

#include <iostream>
#include <vector>

using namespace std;

long long int n, vetor[1000010];


int main()
{
    vetor[4] = 1;
    for(long long int i = 5; i <= 1000001; i++)
        vetor[i] = vetor[i-1] + ((i*(i-1)) - (((i/2)+1)*(i/2)) - ((i+1)*(i-((i/2)+1))));

    while(cin >> n, n >= 3)
    {
        if(n == 3) cout << '0' << endl;
        else cout << vetor[n] << endl;
    }
    return 0;
}
