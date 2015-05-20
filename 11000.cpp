//Anderson Zudio de Moraes
//Victor Crcel Messner
//AC em 12/04/2014

#include <iostream>
#include <math.h>

using namespace std;

int main()
{
    long long int N, macho, bee;
    while(cin >> N, N >= 0)
    {
        if(N == 0)
        {
            cout << "0 1" << endl;
        }
        else if(N == 1)
        {
            cout << "1 2" << endl;
        }
        else
        {
            macho = (1/sqrt(5))*( pow((1 + sqrt(5))/2, N+2) -  pow((1 - sqrt(5))/2, N+2)) - 1;
            bee = (1/sqrt(5))*( pow((1 + sqrt(5))/2, N+3) -  pow((1 - sqrt(5))/2, N+3)) - 1;
            cout << macho << " " << bee << endl;
        }
    }
    return 0;
}
