//Anderson Zudio
//AC em 31/07/2014
#include <iostream>
#include <climits>

using namespace std;

int main()
{
      int N, B, H, W;
      while (cin >> N >> B >> H >> W)
      {
           int r = INT_MAX;
           for (int i = 0; i < H; i++)
           {
                int p;
                cin >> p;
                for (int j = 0; j < W; j++)
                {
                     int a;
                     cin >> a;
                     if (N*p <= B && a >= N && N*p < r)
                     {
                          r = N*p;
                     }
                }
           }
           if (r == INT_MAX)
                cout << "stay home" << endl;
           else
                cout << r << endl;
      }
      return 0;
}
