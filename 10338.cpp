//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 12/04/2014
#include <iostream>
#include<stdio.h>
using namespace std;

int n, p;  long long int C[1001][1001];
long long int repetidos[30],repetidos1[30];

void TrianguloPascal(int n)
{   int i, p;
    C[0][0] = 1; for (p=1; p<=1000; p++) C[0][p]=0;
    for (i=1; i<=n; i++)
    {   C[i][0]=1;
        for (p=1; p<=n; p++) C[i][p]=C[i-1][p]+C[i-1][p-1];
    }
}

int main()
{    //TrianguloPascal(1000);
     long long int val,val1,CONT2=1,casos,aux1,aux,i;
     char a;
     string c;
     cin >> casos;
     while (casos--){
        fill(repetidos,repetidos+28,0);
        fill(repetidos1,repetidos1+28,0);
        val1=1;i=1,val=1;
        cin >> c;
        aux=c.size();
        for(i=0;i<aux;++i){
            aux1=c[i]-'A';
            val1*=i+1;
            if (repetidos[aux1]){
                repetidos1[aux1]++;
                repetidos[aux1]*=repetidos1[aux1];}
            else{
                repetidos[aux1]++;
                repetidos1[aux1]++;}}

        for (i=0;i<27;++i){
            if (repetidos[i]){
                val*=repetidos[i];}}
        cout << "Data set " <<CONT2++ <<": "<< val1/val<<endl;;

        }
return 0;}
