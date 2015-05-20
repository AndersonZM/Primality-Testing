//Anderson Zudio 
//Victor Cracel Messner 
//Julio Cesar
//AC 10/05/2014 

#include <iostream> 
#include <stdio.h> 
#include <string.h> 
using namespace std; 
  
char P[100001], T[100001]; 
int n, m, F[100001], R[256]; 
bool achei; 
  
int KMP(){ 
    int i,j,k; 
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
        if (j == (m-1)){achei = true; return 0;  } 
        if (j == -1) i++; 
        else j = F[j]; 
    } 
    return 0; 
} 
int main() 
{ 
    int t, k; 
    scanf("%d",&t); 
    while(t--) 
    { 
        scanf("%s",T); 
        n = strlen(T); 
        scanf("%d",&k); 
        for(int i = 1; i <= k; i++){ 
            scanf("%s",P); 
            m = strlen(P); 
            achei = false; 
            KMP(); 
            if(achei) printf("y\n"); 
            else printf("n\n"); 
        } 
    } 
    return 0; 
} 