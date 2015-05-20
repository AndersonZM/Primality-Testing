//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014

#include <iostream>

using namespace std;


int main()
{
    long long int N = 1, r;

    while(cin >> N, N >= 0)
    {
        if(N == 0)
        {
            cout << '1' << endl;
            continue;
        }

        r = ((N*(N+1))/2) + 1;
        cout << r << endl;
    }
     return 0;
}
