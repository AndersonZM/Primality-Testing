//Ac em 11/5/2015
//Anderson Zudio
#include <iostream>

using namespace std;

struct node {
    int nf; //Vale 0 se o nó não forma uma palavra caso contrário mostra o número de vezes que está palavra foi inserida na Trie
    node* pont[26]; //Ponteiros para nós de 'a' até 'z'
};

node* gT;
string gS;
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
    int i, gnp, gnpb;
    string c;

    gT = new (node);  gT->nf = 0;  for (i=0; i<26; i++)  gT->pont[i] = NULL; //Cria Trie
    while(cin >> c){
        for(int i = 0; i <= c.size(); i++){
            if(c[i] >= 'A' && c[i] <= 'Z') gS += c[i];
            else if(c[i] >= 'a' && c[i] <= 'z') gS += c[i];
            else{
                if(c.size() == 0) continue;
                InsereTrie();
                gS.clear();
            }
        }
    }
    first = true;
    gS.clear();  ImprimePalavrasTrie(gT, gS);

    return 0;
}


