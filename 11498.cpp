//Anderson Zudio
//AC em 29/07/2014
#include <iostream>

using namespace std;

int main()
{
    int n, i, X, Y, x, y;
    while(cin >> n, n)
    {
        cin >> X >> Y;
        for(i = 0; i < n; i++)
        {
            cin >> x >> y;
            if(x == X) cout << "divisa" << endl;
            else if(y == Y) cout << "divisa" << endl;
            else if(x > X && y > Y) cout << "NE" << endl;
            else if(x > X && y < Y) cout << "SE" << endl;
            else if(x < X && y > Y) cout << "NO" << endl;
            else cout << "SO" << endl;
        }
    }
    return 0;
}
