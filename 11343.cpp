//Anderson Zudio de Moraes
//Julio Cesar
//Victor Cracel Messner
//AC 31/05/2014
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
typedef struct {int x, y;}Ponto;
typedef struct {Ponto p, q;}Segmento;
Ponto pol[100000], pols[100000]; //vetores globais

int SentidoPercurso(Ponto p1, Ponto p2, Ponto p3)
{  long long a, b, c, d;
   a= p2.x - p1.x;  b= p3.y - p1.y;  c=p3.x - p1.x;  d=p2.y - p1.y;
   a = a*b-c*d;
   if (a > 0) return 1;
   else if (a < 0) return -1;
   else return 0;
}
int Intercepta(Segmento s1, Segmento s2)
{  if ((max(s1.p.x,s1.q.x)>=min(s2.p.x,s2.q.x))&&
       (max(s2.p.x,s2.q.x)>=min(s1.p.x,s1.q.x))&&
       (max(s1.p.y,s1.q.y)>=min(s2.p.y,s2.q.y))&&
       (max(s2.p.y,s2.q.y)>=min(s1.p.y,s1.q.y))&&
       (SentidoPercurso(s1.p,s1.q,s2.p)*SentidoPercurso(s1.p,s1.q,s2.q)<=0)&&
       (SentidoPercurso(s2.p,s2.q,s1.p)*SentidoPercurso(s2.p,s2.q,s1.q)<=0))
       return 1;
    else return 0;
}

int main()
{
    int t, m, T;
    bool inter;
    Segmento s[110];
    cin >> t;
    while(t--){
        cin >> m;

        for(int i = 1; i <= m; i++)
            cin >> s[i].p.x >> s[i].p.y >> s[i].q.x >> s[i].q.y;


        T = 0;
        for(int i = 1; i <= m; i++){
            inter = false;
            for(int j = 1; j <= m; j++){
                if(i != j){
                    if(Intercepta(s[i],s[j])) inter = true;
                }
            }
            if(!inter) T++;
        }

        cout << T << endl;
    }
   return 0;
}
