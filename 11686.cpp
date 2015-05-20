//Anderson Zudio de Moraes
//AC em 1/05/2014

#include<iostream>
#include<string.h>
using namespace std;
const int NVM=1000010;

// Declarações para o grafo
typedef struct no* lista;
struct no{ int ve; lista prox;};

int n, m, cpre, cpos, pre[NVM], ot[NVM]; lista A[NVM]; bool ciclo;

//Outras declarações
int i, j, u, w;  lista p;

void Inicializa(int n){
    for (int i=1; i<=n; i++) A[i] = NULL;
    memset(pre, 0, sizeof(pre));
    memset(ot, 0, sizeof(ot));
	cpre = 0;  cpos = n;
}
void OrdTop(int v){
     lista aux = NULL;
     if(pre[v] == 2) ciclo = true;
     else if(pre[v] == 0)
     {
         if(A[v] != NULL) aux = A[v];
         {
            pre[v] = 2;
            while(aux != NULL)
            {
                OrdTop(aux->ve);
                aux = aux->prox;
            }
         }
        ot[cpos--] = v;
        pre[v] = 1;
     }
}

int main(){
    lista aux;
    while(cin >> n >> m, n){
        Inicializa(n);
        for(i = 1;i <= m;i++){
            cin>> u >> w;
			p = new(no);   p->ve = w;   p->prox = A[u];   A[u] = p;
        }
        ciclo = false;
        for (i=1; i<=n; i++)
            if(pre[i] == 0)
                OrdTop(i);
        if(ciclo) cout << "IMPOSSIBLE" << endl;
        else for (i=1; i<=n; i++) cout<<ot[i]<<endl;
    }
    return 0;
}
