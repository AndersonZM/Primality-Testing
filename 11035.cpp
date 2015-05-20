//Anderson Zudio
//AC em 12/05/2015

#include <iostream>
#include <map>

using namespace std;

struct node {
    int nf; //Vale 0 se o nó não forma uma palavra caso contrário mostra o número de vezes que está palavra foi inserida na Trie
    node* pont[52]; //Ponteiros para nós de 'a' até 'z'
};

node* gT;
int gS[100010];
int resp;
map <string, int> okk;

void InsereTrie(int n)
{   int j, d, k, i;  node *p, *q;
    p = gT;
    for(j = n-1; j >= 0; j--)
    {
        d = gS[j]; //converte os caracteres em minúsculo
        if (p->pont[d] == NULL)
        {
            q = new (node);  for (i=0; i<52; i++) q->pont[i] = NULL; resp++;
            q->nf = 0;  p->pont[d] = q;  p = q;
        }
        else p = p->pont[d];
    }
    p->nf++;
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

int Converte(string t){
    char r;

}


int main()
{
    int i, gnp, gnpb, t;
    char naipe[4] = {'C', 'D', 'H', 'S'};
    string value[13] = { "A" , "2" , "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K"};
    string carta;
    int number = 0;
    for(i = 0; i < 4; i++){
        for(int j = 0; j < 13; j++){
            carta = value[j] + naipe[i];
            okk[carta] = number++;
        }
    }


    while(cin >> t, t){
        gT = new (node);  gT->nf = 0;  for (i=0; i<52; i++)  gT->pont[i] = NULL; //Cria Trie
        resp = 0;
        for(i = 0; i < t; i++){
            string palavra;
            cin >> gnp;
            for (int j=0; j<gnp;++j){
                cin >> palavra;
                gS[j] = okk[palavra];
            }
            InsereTrie(gnp);
        }
        cout << resp << endl;
    }

    return 0;
}


