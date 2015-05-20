//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
//AC 16/05/2014

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <sstream>

using namespace std;


int main()
{
    int t, N, M, tamanho[105], aux, aux2, cont = 1, conta_letra1[130], conta_letra2[130];
    string palavras[105], st, linha;
    istringstream ss;
    bool achei, primeiro;
    cin >> t;
    while(t--){
        cin >> N >> M;
        for(int i = 1; i <= N; i++){
            cin >> palavras[i];
            tamanho[i] = palavras[i].size();
        }
        getline(cin, linha);

        //for(int i = 1; i <= N; i++) cout << palavras[i] << " " << tamanho[i] << endl;
        cout << "Case #" << cont++ << ":" << endl;
        for(int i = 1; i <= M; i++){
            ss.clear();
            getline(cin, linha);
            ss.str(linha);
            primeiro = true;
            while(ss >> st){
                aux = st.size();
                for(int j = 1; j <= N; j++){
                    if(tamanho[j] == aux){
                        fill(conta_letra1+95,conta_letra1+130, 0);
                        fill(conta_letra2+95,conta_letra2+130, 0);
                        if(palavras[j][0] == st[0] && palavras[j][aux-1] == st[aux-1]){
                            for(int k = 0; k < aux; k++){
                                aux2 = (int) st[k];
                                conta_letra1[aux2]++;
                            }
                            for(int k = 0; k < aux; k++){
                                aux2 = (int) palavras[j][k];
                                conta_letra2[aux2]++;
                            }
                            achei = true;
                            for(int k = 95; k < 130; k++)
                                if(conta_letra1[k] != conta_letra2[k]){
                                    achei = false;
                                    break;
                                }
                            if(achei){
                                primeiro ? cout << palavras[j] : cout << " " << palavras[j];
                                primeiro = false;
                                break;
                            }
                        }
                    }
                }
            }
            cout << endl;
        }
    }
    return 0;
}
