//Anderson Zudio
//AC em 29/07/2014
#include <iostream>

using namespace std;

int main()
{
    int t, x, y, a, cont = 1;
    while(cin >> t, t){
        x = y = 0;
        for(int i = 1; i <= t; i++)
        {
            cin >> a;
            a ? x++ : y++;

        }
        cout << "Case " << cont++ << ": " << x - y << endl;
    }
    return 0;
}
