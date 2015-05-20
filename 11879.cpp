//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014
#include <iostream>
#include <string>
using namespace std;
int a, b, n, K;

string r;

int main()
{
    int i;
    r = "1";
    while(cin >> r, r[0] != '0')
    {

        a = 0;
        i = 0;
        while(r[i])
        {

            a *= 10;
            a += r[i] - '0';
            a = a % 17;
            i++;
        }




        if(a == 0) cout << 1 << endl;
        else cout << 0 << endl;

    }
    return 0;
}
