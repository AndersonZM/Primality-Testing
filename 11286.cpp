//Anderson Zudio
//AC em 12/05/2015

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

struct node {
    int nf; //Vale 0 se o nó não forma uma palavra caso contrário mostra o número de vezes que está palavra foi inserida na Trie
    node* pont[401]; //Ponteiros para nós de 'a' até 'z'
};

node* gT;
node* resp2[10001];
vector<int> gS;
int resp;
map <string, int> okk;

void InsereTrie(int n)
{   int j, d, k, i;  node *p, *q;
    p = gT;
    for(j = 0; j < 5; j++)
    {
        d = gS[j]; //converte os caracteres em minúsculo
        if (p->pont[d] == NULL)
        {
            q = new (node);  for (i=0; i<401; i++) q->pont[i] = NULL;
            q->nf = 0;  p->pont[d] = q;  p = q;
        }
        else p = p->pont[d];
    }
    p->nf++;
    resp2[n] = p;
    resp = max(resp, p->nf);

}

void ImprimePalavrasTrie(node* p, string t)
{    int d;
     if (p->nf)  cout << t << " (" << p->nf << ")" << endl; //Verifica se no node tem uma palavra, se sim imprime ela e o numero de ocorrencia dela
     for (d=0; d < 26; d++)
        if (p->pont[d] != NULL)
        {
            t += 'a' + d;
            ImprimePalavrasTrie(p->pont[d], t);
            t.erase(t.length() - 1); //t.pop_back() //c++11
        }
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
    int i, gnp, gnpb, t;

    while(cin >> t, t){
        gT = new (node);  gT->nf = 0;  for (i=0; i<401; i++)  gT->pont[i] = NULL; //Cria Trie
        resp = -1;
        for(i = 0; i < t; i++){
            for (int j=0; j < 5;++j){
                int y;
                cin >> y;
                gS.push_back(y-100);
            }
            sort(gS.begin(), gS.end());
            InsereTrie(i);
            gS.clear();
        }
        int aux = 0;
        for(i = 0; i < t; i++){
            if(resp2[i]->nf == resp) aux++;
        }

        cout << aux << endl;
    }

    return 0;
}


