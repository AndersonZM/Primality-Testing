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
int tam, n;
bool first;
int resp;

void InsereTrie(void)
{   int j, d, k, i;  node *p, *q;
    p = gT;
    for(j = 0; j < gS.length(); j++)
    {
        d = tolower(gS[j]) - 'a'; //converte os caracteres em minúsculo
        if (p->pont[d] == NULL)
        {
            q = new (node);  for (i=0; i<26; i++) q->pont[i] = NULL;
            q->nf = 1;  p->pont[d] = q;  p = q;
        }
        else{
            p = p->pont[d];
            p->nf++;
        }
    }
}

void ImprimePalavrasTrie(node* p, string t, int ant)
{    int d;
     if (p->nf != ant){ resp += p->nf;} //Verifica se no node tem uma palavra, se sim imprime ela e o numero de ocorrencia dela
     for (d=0; d < 26; d++)
        if (p->pont[d] != NULL)
        {
            t += 'a' + d;
            ImprimePalavrasTrie(p->pont[d], t, p->nf);
            t.pop_back(); //c++11
        }
}


int main()
{
    int i, gnp, gnpb, t = 0;
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);

    while(cin >> n){
        gT = new (node);  gT->nf = 0;  for (i=0; i<26; i++)  gT->pont[i] = NULL; //Cria Trie
        for(int i = 0; i < n; i++){
                cin >> gS;
                InsereTrie();
        }
        resp = 0;
        gS.clear();
        ImprimePalavrasTrie(gT, gS, 0);
        double r = resp / (double) n;
        cout << r << endl;
    }

    return 0;
}


