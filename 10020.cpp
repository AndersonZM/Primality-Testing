//Anderson Zudio
//Bruno Masquio
//AC em 1/7/2014
#include <iostream>
#include <algorithm>
#include <vector>
#define inf 2147483646

using namespace std;

vector <pair<int, int> > vetor;

bool cmp(pair<int, int> a, pair<int, int> b){ return a.first < b.first; }

void resolve(int M){
    vector <pair <int, int> > R;
    auto q = 0;
    auto p = 0;
    vetor.push_back({inf, inf});
    for(auto i = 1; i < vetor.size(); i++){
        if(vetor[i].first > p){
            R.push_back({vetor[q].first, vetor[q].second});
            p = vetor[q].second; q = i;
            if(p >= M) break;
        }
        else if(vetor[i].second > vetor[q].second) q = i;
    }

    cout << R.size() << endl;
    for(auto i = R.begin(); i != R.end(); i++){
        cout << i->first << " " << i->second;
        if(i != R.end()-1) cout << endl;
    }


}

int main(){
    int t, M, c, f;
    cin >> t;
    while(t--){
        cin >> M;
        vetor.clear();
        vetor.push_back({-inf, -inf});
        while(cin >> c >> f, c != 0 || f != 0)
            vetor.push_back({c, f});

        sort(vetor.begin(), vetor.end(), cmp);
        auto cobre = true;
        auto p = 0;

        for(auto i = vetor.begin()+1; i != vetor.end(); i++){
            if(p < M){
                if(i->first > p) cobre = false;
                else if(i->second > p) p = i->second;
            }
        }
        if(p < M) cobre = false;
        if(cobre) resolve(M);
        else cout << '0';
        if(t) cout << endl << endl;
    }
    return 0;
}
