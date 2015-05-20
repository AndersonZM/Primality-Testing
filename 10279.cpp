//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 24/03/2014

#include <iostream>
#include <stdio.h>

using namespace std;

int main(){
    int testes,n,i,j,booleano,r[12][12];
     char a[12][12],b[12][12];
    cin >> testes;
    while(testes--){
        cin >> n;
        booleano=0;

        for (i=0;i<12;++i){
            for (j=0;j<12;++j){
                r[i][j]=0;}}

        for(i=1;i<=n;i++)
            for(j=1;j<=n;j++)
                cin >> a[i][j];

        for(i=1;i<=n;i++)
            for(j=1;j<=n;j++)
            {
                cin >> b[i][j];
                if(a[i][j] == '*' && b[i][j] == 'x')
                    booleano = 1;
            }


        for(i=1;i<=n;i++)
            for(j=1;j<=n;++j){
                if(a[i][j]=='*'){
                    r[i-1][j-1]+=1;
                    r[i-1][j]+=1;
                    r[i-1][j+1]+=1;
                    r[i][j-1]+=1;
                    r[i][j+1]+=1;
                    r[i+1][j-1]+=1;
                    r[i+1][j]+=1;
                    r[i+1][j+1]+=1;}}

    for (i=1;i<=n;++i){
        for (j=1;j<=n;++j){
            if(booleano && a[i][j] == '*')
                cout << '*';
            else if(b[i][j] == 'x')
                cout << r[i][j];
            else
                cout << '.';

        }
        cout<<endl;}
    if (testes){
        cout<<endl;}}
return 0;
}
