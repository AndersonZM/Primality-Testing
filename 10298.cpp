//Ac em 10/05/2014
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;

string P, T;
int n, m, F[1000001], R[256], cont;
bool achei;

void KMP(){
    int i,j,k, elegante, cont;
    cont = elegante = 0;
    //Preprocessamento
    F[0] = -1;  j = -1;  i = 0;
    while (i < (m-1)) {
        if (P[i+1] == P[j+1]){ j++; F[++i] = j; }
        else if (j = -1)  F[++i] = -1;
        else j = F[j];
    }
    //Busca
    i=-1;  j=-1;
    while ((i-j) <= (n-m)){
        while (j < (m-1)){
            if (P[j+1] == T[i+1]){ i++;  j++;}
            else break;
        }
        if (j == (m-1)) {
            if((i+1 - elegante) == m){
                    elegante = i+1;
                    cont++;
                    if(elegante == n){
                        achei = true;
                        cout << cont << endl;
                    }
            }
            else break;
        }
        if (j == -1) i++;
        else j = F[j];
    }
}

int main()
{
    int r, aux;
    while(cin >> T, T[0] != '.')
    {
        n = T.size();
        P = "";
        achei = false;
        m = 1;
        for(int l = 0; l < n/2; l++){
            P += T[l];
            if((n % m) == 0){
                KMP();
            }
            m++;
            if(achei) break;
        }
        if(!achei) cout << "1" << endl;
    }
    return 0;
}

