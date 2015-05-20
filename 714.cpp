//Anderson Zudio
//AC 29/07/2014
#include<iostream>
#include<string>
#include <vector>

using namespace std;

int m, k;
long long s, minimo;

vector <long long int> numeros;
vector <bool> vis(510);


int divide(long long int M){
    fill(vis.begin(), vis.begin()+510, false);
    int cont = 0;
    int pos = m-1;
    while(pos >= 0){
        long long int soma = 0;
        auto aux = true;
        while(pos >= 0 && soma + numeros[pos] <= M){
            aux = false;
            soma += numeros[pos];
            pos--;
        }
        if(aux) return k + 1;
        if(pos >= 0) vis[pos] = true;
        cont++;
    }
    return cont++;
}

long long int resolve(){
    long long int l = minimo, r = s, meio;
    while(l < r){
        meio = (l + r) >> 1; //10^(l+r)
        if( divide(meio) <= k)
            r = meio;
        else
            l = meio + 1;
    }
    return r;
}

inline void output(){

}

int main(){
    int t;
    cin >> t;
    while(t--){
        cin >> m >> k;
        s = minimo = 0;
        numeros.clear();
        long long numero;
        for(auto i = 0; i < m; i++){
            cin >> numero;
            numeros.push_back(numero);
            s += numero;
            if(numero > minimo) minimo = numero;
        }

        auto resp = resolve();

        auto cont = divide(resp);
        for(int i = 0; i < m-1 && cont < k; i++)
            if(!vis[i]){
                vis[i] = true;
                cont++;
            }
        for(int i = 0; i < m; i++){
            if(i) cout << " " << numeros[i];
            else cout << numeros[i];
            if(vis[i]){
                cout << " /";
            }
        }
        cout << endl;
    }
    return 0;
}
