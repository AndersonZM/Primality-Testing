//Ac em 10/05/2014 
//Anderson Zudio 
//Victor Cracel Messner 
  
#include <iostream> 
#include <stdio.h> 
#include <string.h> 
#include <string> 
#include <algorithm> 
using namespace std; 
  
//char P[1000001], T[1000001]; 
string k[25], e[25], P, T; 
int n, m, F[1000001], r[256]; 
  
int KMP(){ 
    int i,j,k, cont = 0; 
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
        if (j == (m-1)){ 
           if(i >= n){ 
                if((!(T[i-m-1] >= 'A' && T[i-m-1] <= 'Z')) && (!(T[i+1] >= 'A' && T[i+1] <= 'Z'))){ 
                    cont++; 
                } 
            } 
            else{ 
                if(!(T[i+1] >= 'A' && T[i+1] <= 'Z')){ 
                    cont++; 
  
                } 
            }//cont++; cout << "debug " << i << " " << m <<" " <<T[i] << endl; 
        } 
        if (j == -1) i++; 
        else j = F[j]; 
    } 
  
    return cont; 
} 
  
int main() 
{ 
    int i, j, K, E, melhor, cont = 0; 
    string lixo; 
    while(scanf("%d %d", &K, &E) != EOF){ 
            getline(cin,lixo); 
        for(i = 0; i < K; i++){ 
            getline(cin ,k[i]); 
        } 
        melhor = 0; 
        for(i = 0; i < E; i++){ 
            getline(cin ,e[i]); 
            T = e[i]; 
            n = T.size(); 
            transform(T.begin(),T.end(),T.begin(),::toupper); 
            r[i] = 0; 
            for(j = 0; j < K; j++){ 
                    P = k[j]; 
                    m = P.size(); 
                    transform(P.begin(),P.end(),P.begin(),::toupper); 
                    r[i] += KMP(); 
            } 
            melhor = max(melhor, r[i]); 
        } 
        cout << "Excuse Set #" << ++cont << endl; 
        for(i = 0; i < E; i++){ 
            if(r[i] == melhor) cout << e[i] << endl; 
        } 
        cout << endl; 
    } 
    return 0; 
} 