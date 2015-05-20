//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//AC em 21/05/2014

#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;

string P, matriz[103];
int n, m, F[1000001], R[256];
bool achei;
pair <int, int> cini;
pair <int, int> cfim;

int main()
{
    int n, m, i, j, k, l, w;
    string lixo;
    cin >> n;
    getline(cin, lixo);
    for(int i = 0; i < n; i++){
        getline(cin, matriz[i]);
    }

    while(cin >> P, P[0] != '0'){
            m = P.size();
            achei = false;


            for(i = 0; i < n; i++){
                for(j = 0; j < n; j++){
                    if(matriz[i][j] == P[0]){
                        //Processamento começa na direita e roda no sentido horário
                        if(j+m-1 < n){
                            k = j+1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[i][k++] != P[l++]){
                                     achei = false;
                                     break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i+1,j+m);
                            }
                        }

                        if(!achei && i+m-1 < n && j+m-1 < n){
                            k = j+1;
                            w = i+1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w++][k++] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i+m,j+m);
                            }
                        }

                        if(!achei && i+m-1 < n){
                            w = i+1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w++][j] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i+m,j+1);
                            }
                        }

                        if(!achei && i+m-1 < n && j-m+1 >= 0){
                            k = j-1;
                            w = i+1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w++][k--] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i+m,j-m+2);
                            }
                        }

                        if(!achei && j-m+1 >= 0){
                            k = j-1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[i][k--] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i+1,j-m+2);
                            }
                        }

                        if(!achei && i-m+1 >= 0 && j-m+1 >= 0){
                            k = j-1;
                            w = i-1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w--][k--] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i-m+2,j-m+2);
                            }
                        }

                        if(!achei && i-m+1 >= 0){
                            w = i-1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w--][j] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i-m+2,j+1);
                            }
                        }

                        if(!achei && i-m+1 >= 0 && j+m-1 < n){
                            k = j+1;
                            w = i-1;
                            l = 1;
                            achei = true;
                            while(l < m){
                                if(matriz[w--][k++] != P[l++]){
                                    achei = false;
                                    break;
                                }
                            }
                            if(achei){
                                cini = pair<int, int>(i+1,j+1);
                                cfim = pair<int, int>(i-m+2,j+m);
                            }
                        }

                    }
                    if(achei) break;
                }
                if(achei) break;
            }

            if(achei) cout << cini.first << "," << cini.second << " " << cfim.first << "," << cfim.second << endl;
            else cout << "Not found" << endl;


    }


    return 0;
}
