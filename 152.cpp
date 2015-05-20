//Anderson Zudio
//Victor Cracel
//AC em 27/08/2014

#include<stdio.h>
#include <iostream>
#include <cmath>

using namespace std;

int main(){
    long long int dist[10], x[5010],y[5010],z[5010];
    double modulo[5010];
    long long int a,b,c;
    double sub, sub2;
    int caso=0,i,j;
    for (i=0;i<10;++i){
        dist[i]=0;}
    cin >> a >> b >> c;
    while (a||b||c){
        ++caso;
        x[caso]=a;y[caso]=b;z[caso]=c;
        cin >> a >> b >> c;}
    for (i=1;i<=caso;++i){
    	sub=1<<31-1;
        for (j=1;j<=caso;++j){
        	if (i==j) continue;
            sub2 = sqrt(((x[i]-x[j])*(x[i]-x[j])) + ((y[i]-y[j])*(y[i]-y[j])) + ((z[i]-z[j])*(z[i]-z[j])));
            if (sub > sub2) sub = sub2;
            }
        if(sub < 10)dist[(int)sub]++;
            }
    for (i=0;i<10;++i){
    	printf ("%4d",dist[i]);
    }
    printf("\n");
    



    return 0;
}