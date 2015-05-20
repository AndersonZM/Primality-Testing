//Anderson Zudio
//Bruno Masquio
//AC 27/05/2014

#include<iostream>
#include <string.h>
#define MAX 2147483485
using namespace std;
const int NVM=150;

/* A Solução deste problema é feita com o algoritmo de Dijkstra, a idéia é minimizar os caminhos
para todos os vértices do grafo em relação ao vértice de saída. Para isso
o programa recebe todas as arestas invertidas e faz a solução com Dijkstra para o vértice
que o problema chama de saída. Após obter a resposta em Di, o program varre Di e ve quem o
vertice de saída pode alcançar dado T (o peso em unidade de tempo)
*/


// Declarações para o grafo
int n, m, E[NVM][NVM], Di[NVM], pre[NVM];
//Outras declarações
int i, j, u, w, p, z;

void Inicializa(){
    memset(E, -1, sizeof(E));
    memset(pre, 0, sizeof(pre));
    fill(Di, Di+NVM, MAX);
}

int Solucao(int e, int temp){
	int cont;
	//Dijkstra
	Di[e] = 0;
	
	for(i = 1; i <= n; i++){
		for(j = 1; j <= n; j++){
			if(!pre[j]){
				w = j;
				break;
			}
		}
		
		for(j = w+1; j <= n; j++){
			if((!pre[j]) && (Di[j] < Di[w])){
				w = j;
			}
		}
		
		pre[w] = 1;
		
		for(z = 1; z <= n; z++){
			if((pre[z] == 0) && (E[w][z] > -1)){
				if(Di[w] + E[w][z] < Di[z]) Di[z] = Di[w] + E[w][z];
			}
		}
	}
	cont = 0;
	//A contagem de quantos vértices a saída alcança
	for(i = 1; i <= n; i++)
		if(Di[i] <= temp) cont++;
		
	return cont;
}

int main(){
	int t, exi, temp;
	cin >> t;

	while(t--){
		cin >> n >> exi >> temp >> m;
	    Inicializa();
	    
	    for(i=1; i<=m; i++){
		    cin >> u >> w >> p;
		//Recebendo as arestas invertidas
		    E[w][u] = p;
	    }
		//O retorno da solução é quantos vértices o de saída alcança dado Temp
	    cout << Solucao(exi, temp) << endl;
	    if(t) cout << endl;
    }
	return 0;
}