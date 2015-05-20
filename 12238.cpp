//Anderson Zudio
//Raquel Marcolino
//Ac em 11/04/2015

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cmath>
#define NMAX 500000
#define addFilho filhos.push_back
#define nFilhos filhos.size()
#define LCA(A,B) E[RMQ(H[A],H[B])]
using namespace std;

//ATENÃÃO: 	LCA nÃ£o Ã© uma funÃ§Ã£o, as chamadas sÃ£o diretamente convertidas
//			em RMQ na macro acima.

//Estrutura da Ã¡rvore: array de nÃ³s com Ã­ndices dos filhos
struct no {
	vector<long long int> filhos;
};

no T[NMAX];
//Estruturas e variÃ¡veis auxiliares para a DFS e o AE
long long int E[3*NMAX], L[3*NMAX], H[NMAX], M[3*NMAX][50], J[NMAX], ii;

//ConstruÃ§Ã£o do circuito euleriano, e dos arrays L e H
void dfs(int n, int nv) {
	int f;
	E[ii]=n; L[ii]=nv;
	H[n]=ii++; //Primeira ocorrÃªncia do nÃ³ n
	for(f=0;f<T[n].nFilhos;++f) { //Visita todos os filhos de n
		dfs(T[n].filhos[f], nv+1);
		E[ii]=n; L[ii++]=nv;
	}
}

//RMQ com Array Esparso
int menor(int a, int b) {
	return (L[b]>L[a]) ? a : b;
}

void arrayEsparso(int N) {
	N *= 2;
	int i, j;
    for (i = 0; i < N; ++i) M[i][0] = i;
    for (j = 1; 1<<j <= N; ++j)
        for (i = 0; i+(1<<j)-1 < N; ++i) {
        	M[i][j] = menor(M[i][j-1],M[i+(1<<(j-1))][j-1]);
        }
}

int RMQ(int c, int f) {
	if (c>f) {int aux = c; c=f; f=aux; } //Troca para evitar consultas invÃ¡lidas
	int tam = f-c+1;
	int k;
	for (k=1; 1 << k <= tam; ++k); --k; //cÃ¡lculo de piso(log(tam))
	return menor(M[c][k], M[f-(1<<k)+1][k]);
}

////////////////////////////////
int main() {
	int girls,who,say,N,i;
	int att,plx;
	int mono,b,lock;
	//Attention, please!
	while(scanf(" %d",&N)){
		if(!N) break;
		memset(H,0,sizeof(H));
		J[1]=0;
		for(i=1;i<N;++i){
			scanf(" %d %d",&girls,&say);
			T[girls+1].addFilho(i+1);
			J[i+1]=J[girls+1]+say;
		}
		ii=0; dfs(1,0); arrayEsparso(N);
		scanf(" %d",&who);
		lock=1;
		for(i=0;i<who;++i){
			scanf(" %d %d",&att,&plx);
			++att;
			++plx;
			(lock)?lock=0:printf(" ");
			printf("%lld",J[att]+J[plx]-2*J[LCA(att,plx)]);
		}
		printf("\n");
		for(i=0;i<N;++i){
			T[i].filhos.clear();
		}
	}
	return 0;
}
