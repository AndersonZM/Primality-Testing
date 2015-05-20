//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014

#include <iostream>
#include <math.h>

using namespace std;

int main()
{
    long long int D = 1, r;
    int cont = 1;

    while(cin >> D, D > 0)
    {
        r = ceil((3 + (sqrt(9 + (8*D)))) / 2);
        cout << "Case " << cont++ << ": " <<r << endl;
    }
    return 0;
}
