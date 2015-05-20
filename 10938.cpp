//Anderson Zudio
//AC em 25/04/2015

#include<iostream>
#include<vector>
#include<string.h>
#define MAX 5010

using namespace std;

vector <int> grafo[MAX];
bool vis[MAX];
int c[MAX], dist, x;

void dfs(int u, int k){
    vis[u] = true;
    c[k] = u;

    if(u == x){
        dist = k;
        return;
    }

    for(auto i = grafo[u].begin(); i != grafo[u].end(); i++){
        if(vis[*i])
            continue;
        dfs(*i, k+1);
        if(dist != -1)
            return;
    }
}

int main(){
    int n, u, v;
    while(cin >> n, n){

        for(int i = 0; i <= n; i++) grafo[i].clear();

        for(int i = 0; i < n-1; i++){
            cin >> u >> v;
            grafo[u].push_back(v);
            grafo[v].push_back(u);
        }

        int l;
        cin >> l;

        for(int i = 0; i < l; i++){
            memset(vis, 0, sizeof(vis));
            dist = -1;

            cin >> u >> x;

            dfs(u, 0);
            if(dist % 2 == 0 )  cout << "The fleas meet at " << c[dist/2] << "." << endl;
            else{
                auto p = c[dist/2];
                auto r = c[dist/2+1];

                cout << "The fleas jump forever between ";
                (p > r) ? cout << r << " and " << p << "." : cout << p << " and " << r << ".";
                cout << endl;
            }
        }
    }
}
