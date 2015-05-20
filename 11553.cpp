//Anderson Zudio
//Bruno Masquio
//AC em 1/7/2014
#include <algorithm>
#include <iostream>
#include <vector>
#define INF 2147483646

using namespace std;

vector <int> meuvector;
int M[9][9];

int main(){
    int t, n;
    cin >> t;
    for(auto i = 0; i < 8; i++)
        meuvector.push_back(i);
    while(t--){
        cin >> n;
        for(auto i = 0; i < n; i++)
            for(auto j = 0; j < n; j++)
                cin >> M[i][j];

        auto resp = INF;
        do{
            auto soma = 0;
            for(auto i = 0, j = 0; i < n; i++, j++)
                soma += M[i][meuvector[j]];
            resp = min(soma, resp);
        }while(next_permutation(meuvector.begin(), meuvector.begin()+n));

        cout << resp << endl;
    }

    return 0;
}
