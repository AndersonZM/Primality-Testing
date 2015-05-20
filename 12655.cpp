//Anderson Zudio
//30/04/2015
#include<iostream>
#include<algorithm>
#include <string.h>
#define nFilhos filhos.size()
#define LCA(A,B) E[RMQ(H[A],H[B])]
using namespace std;
const int NVM2=100010;
const int INF=1000000000;
struct aresta {
    int u; int v; int p;
};

struct no {
	vector< pair<int, int> > filhos;
};

no T[NVM2];

int E[NVM2], L[NVM2], H[NVM2], M[NVM2][50], P[NVM2][50], ii, O[NVM2], pre[NVM2];

//Construção do circuito euleriano, e dos arrays L e H
void dfs(int n) {
	pre[n] = 1;
	for(auto f = T[n].filhos.begin(); f != T[n].filhos.end(); ++f) { //Visita todos os filhos de n
        if(!pre[f->first]){
            M[f->first][0] = n;
            P[f->first][0] = f->second;
            L[f->first] = L[n]+1;
            O[f->first] = f->second;
            dfs(f->first);
        }
	}
}

//RMQ com Array Esparso
int menor(int a, int b) {
	return (L[b]>L[a]) ? a : b;
}

void arrayEsparso(int N) {
	int i, j;
    /*for (i = 0; i <= N; ++i){
        M[i][0] = i;
        P[i][0] = O[i];
    }*/
    for (j = 1; (1<<j) < N; ++j){
        for (i = 1; i <= N; ++i) {
            if(M[i][j-1] != -1){
                M[i][j] = M[M[i][j-1]][j-1];
                P[i][j] = min(P[i][j-1], P[M[i][j-1]][j-1]);
            }
        }
    }
}

int RMQ(int c, int f) {
    int resp = INF;
	if (L[c] < L[f]) {int aux = c; c=f; f=aux; } //Troca para evitar consultas inválidas
	int tam = f-c+1;
	int k;
	for (k=0; (1 << k) <= L[c]; ++k); --k; //cálculo de piso(log(tam))


	if(L[c] != L[f]){
        for(int i = k; i >= 0; i--){
            if(L[c] - (1 << i) >= L[f]){
                    resp = min(resp, P[c][i]);
                    c = M[c][i];
            }
        }
	}

	if(c == f) return resp;

	for(int i = k; i >= 0; i--){
        if(M[c][i] != -1 && M[c][i] != M[f][i]){
            resp = min(min(P[c][i], resp), P[f][i] );
            c = M[c][i];
            f = M[f][i];
        }
	}

	resp = min(resp, min(O[c], O[f]));

	return resp;
}

// Declarações para o grafo
int n, m, nt, ct, Pai[NVM2];
aresta EE[NVM2];

int compara(aresta e, aresta f){
    return e.p > f.p;
}
int UFFind(int p){
    if (Pai[p] != p)  Pai[p] = UFFind(Pai[p]);
    return Pai[p];
}
void UFUnion(int p, int q){
    if (Pai[p] > Pai[q]) Pai[Pai[p]] = Pai[q];
    else Pai[Pai[q]] = Pai[p];
}

void Kruskal(){
   int i;
   sort(EE+1, EE+m+1, compara);   nt=0;  ct = 0;
   for(i = 1; i <= n; i++){
       Pai[i] = i;
   }
   for(i = 1; i <= m; i++){
       UFFind(EE[i].u);  UFFind(EE[i].v);
       if(Pai[EE[i].u] != Pai[EE[i].v]){
	       UFUnion(EE[i].u, EE[i].v);
                nt++;
                T[EE[i].u].filhos.push_back( {EE[i].v, EE[i].p} ); ct = ct + EE[i].p;
                T[EE[i].v].filhos.push_back( {EE[i].u, EE[i].p} );
       }
   }
}

int main(){
    int u, v, p, s;
	while (cin >> n >> m >> s){
	    // cout<<endl<<"Grafo com n m = ";  cin >>n>>m;
	    for(int i = 0; i <= n; i++) T[i].filhos.clear();
	    memset(EE, 0, sizeof(EE));
	    memset(M, -1, sizeof(M));
        memset(P, -1, sizeof(P));
        memset(L, -1, sizeof(L));
        memset(pre, 0, sizeof(pre));
	    //cout<<"Arestas e pesos:"<<endl;
	    for(int i=1; i<=m; i++){
		    cin >>EE[i].u>>EE[i].v>>EE[i].p;
	    }
	    Kruskal();
        //cout<<"AGM:"<<endl;
        /*for(int i = 1; i <= n; i++)
        	cout << EE[i].u << " " << EE[i].v << " " << EE[i].p << endl;
        cout << endl;
        cout<<"custo minimo: "<<ct<<endl;*/
        L[1] = 0;
        dfs(1); arrayEsparso(n);

        //for(int i = 1; i <= n; i++)
        //    cout << O[i] << endl;

        /*for(int i = 1; i <= n; i++){
            for(auto j = T[i].filhos.begin(); j != T[i].filhos.end(); j++){
                cout << i << " " << j->first << " " << j->second << endl;
            }
        }*/
        for(int i = 1; i <= s; i++){
            cin >> u >> v;
            cout << RMQ(u , v) << endl;
        }

	}
	return 0;
}


