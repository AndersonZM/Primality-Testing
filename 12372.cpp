//Anderson Zudio
//AC em 29/07/2014
#include <iostream>

using namespace std;

int main()
{
    int t, l, h, w;
    cin >> t;
    for(int i = 1; i <= t; i++)
    {
        cin >> l >> w >> h;
        (l <= 20 && w <= 20 && h <= 20) ? cout << "Case " << i << ": good" << endl : cout << "Case " << i << ": bad" << endl;
    }
    return 0;
}
