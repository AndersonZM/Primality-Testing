//Ac em 14/05/2015
//Anderson Zudio
#include<stdio.h>
#include<iostream>
#include<cstring>
#include<string>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

struct node {
    int nf; //Vale 0 se o nó não forma uma palavra caso contrário mostra o número de vezes que está palavra foi inserida na Trie
    node* pont[26]; //Ponteiros para nós de 'a' até 'z'
};

node* gT;
string gS, c;
int tam, cortes[300005];
bool first;

void InsereTrie(void)
{   int j, d, k, i;  node *p, *q;
    p = gT;
    for(j = 0; j < gS.length(); j++)
    {
        d = tolower(gS[j]) - 'a'; //converte os caracteres em minúsculo
        if (p->pont[d] == NULL)
        {
            q = new (node);  for (i=0; i<26; i++) q->pont[i] = NULL;
            q->nf = 0;  p->pont[d] = q;  p = q;
        }
        else p = p->pont[d];
    }
    p->nf++;
}

void ImprimePalavrasTrie(node* p, string t)
{    int d;
     if (p->nf){ if(!first) cout << t << endl; else first = false;} //Verifica se no node tem uma palavra, se sim imprime ela e o numero de ocorrencia dela
     for (d=0; d < 26; d++)
        if (p->pont[d] != NULL)
        {
            t += 'a' + d;
            ImprimePalavrasTrie(p->pont[d], t);
            t.pop_back(); //c++11
        }
}

int corta(int k) {
    auto p = gT;
    if(k == tam) return 1;

    if(cortes[k] != -1 ) return cortes[k];
    cortes[k] = 0;

    for(int j = k; j < tam; j++){
        p = p->pont[ c[ j ] - 'a' ];
        if(p == NULL )break;
        if( p->nf ) cortes[k] = ( cortes[k] + corta(j+1) ) % 20071027;
     }


     return cortes[k] % 20071027 ;
}


bool BuscaPalavraTrie(node *p, string t){
    int i, k, pl;
    k = t.size();
    for(i = 0; i < k; i++){
        pl = t[i] - 'a';
        if(p->pont[pl] != NULL){
            p = p->pont[pl];
        }
        else return false;
    }

    if(p->nf) return true;
    else return false;
}


int main()
{
    int i, gnp, gnpb, t = 0;


    while(cin >> c){
        gT = new (node);  gT->nf = 0;  for (i=0; i<26; i++)  gT->pont[i] = NULL; //Cria Trie
        cin >> gnp;
        for(int i = 0; i < gnp; i++){
                cin >> gS;
                InsereTrie();
        }
        tam = c.size();
        memset(cortes, -1, sizeof cortes);
        cout << "Case " << ++t << ": " << corta(0) << endl;
    }

    return 0;
}


