//Anderson Zudio
//AC em 3/7/2014

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main()
{
    int n;
    cin >> n;
    for(int i = 1; i <= n; i++)
    {
        vector <int> lista(3, 0);
        cin >> lista[0] >> lista[1] >> lista[2];
        sort(lista.begin(), lista.end());
        cout << "Case " << i << ": " << lista[1] << endl;
    }
    return 0;
}
