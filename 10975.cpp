//Ac em 12/5/2015
//Anderson Zudio
#include <iostream>
#include <map>

using namespace std;

struct node {
    int nf; //Vale 0 se o nó não forma uma palavra caso contrário mostra o número de vezes que está palavra foi inserida na Trie
    node* pont[26]; //Ponteiros para nós de 'a' até 'z'
};

node* gT;
string gS;
bool first;
map <string, int> resp;

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
    p->nf = 1;
}

void ImprimePalavrasTrie(node* p, string t)
{    int d;
     if (p->nf > 1){
            cout << t << " " << p->nf-1 << endl;
            p->nf = 1;
     } //Verifica se no node tem uma palavra, se sim imprime ela e o numero de ocorrencia dela
     for (d=0; d < 26; d++)
        if (p->pont[d] != NULL)
        {
            t += 'a' + d;
            ImprimePalavrasTrie(p->pont[d], t);
            t.pop_back(); //c++11
        }
}

int BuscaPalavraTrie(node *p, string t){
    int i, k, pl;
    bool first = true;
    k = t.size();
    for(i = 0; i < k; i++){
        pl = t[i] - 'a';
        if(p->pont[pl] != NULL){
            p = p->pont[pl];
            first = false;
        }
        else{
                if(first) return -1;
                else return 0;
        }
    }

    if(p->nf){  resp[t] = p->nf; p->nf++; return 1;}
    else return 0;
}
int main()
{
    int i, gnp, gnpb, t, cont1 = 1, cont2 = 1;
    string c[103];
    cin >> t;
    while(t--){
        gT = new (node);  gT->nf = 0;  for (i=0; i<26; i++)  gT->pont[i] = NULL; //Cria Trie
        int i, m, n;
        cin >> i;
        for(; i > 0; i--){
            cin >> gS;
            if(gS.size() <= 100) InsereTrie();
        }
        cout << "Test Case #" << cont1++ << endl;
        cont2 = 1;
        cin >> i;
        for(; i > 0; i--){
            cin >> n >> m;
            for(int u = 0; u < n; u++){
                    cin >> c[u];
            }
            int aux;
            for(int u = 0; u < n; u++){
                for(int v = 0; v < m; v++){

                            node *p;
                            int pl;

                            //esquerda
                            gS.clear();
                            p = gT;
                            for(int k = v; k < m; k++){
                                pl = c[u][k] - 'a';
                                gS += c[u][k];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //direita
                             gS.clear();
                             p = gT;
                            for(int k = v; k >= 0; k--){
                                pl = c[u][k] - 'a';
                                gS += c[u][k];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //baixo
                             gS.clear();
                             p = gT;
                            for(int k = u; k < n; k++){
                                pl = c[k][v] - 'a';
                                gS += c[k][v];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //cima
                            gS.clear();
                            p = gT;
                            for(int k = u; k >= 0; k--){
                                pl = c[k][v] - 'a';
                                gS += c[k][v];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //diagonal esquerda baixo
                            gS.clear();
                            p = gT;
                            for(int k = u, l = v; k < n && l < m; k++, l++){
                                pl = c[k][l] - 'a';
                                gS += c[k][l];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //diagonal esquerda cima
                            p = gT;
                            gS.clear();
                            for(int k = u, l = v; k < n && l >= 0; k++, l--){
                                pl = c[k][l] - 'a';
                                gS += c[k][l];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                           //diagonal direita baixo
                            p = gT;
                            gS.clear();
                            for(int k = u, l = v; k >= 0 && l < m; k--, l++){
                                pl = c[k][l] - 'a';
                                gS += c[k][l];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }

                            //diagonal direita cima
                            p = gT;
                            gS.clear();
                            for(int k = u, l = v; k >= 0 && l >= 0; k--, l--){
                                pl = c[k][l] - 'a';
                                gS += c[k][l];
                                if(p->pont[pl] != NULL){
                                    p = p->pont[pl];
                                    if(p->nf){
                                        resp[gS]++;
                                    }
                                }
                                else break;
                            }
                    }
                }

            cout << "Query #" << cont2++ << endl;
            gS.clear();
            for(auto q = resp.begin(); q != resp.end(); q++){
                cout << q->first << " " << q->second << endl;
            }
            resp.clear();
        }
        cout << endl;
    }

    return 0;
}


