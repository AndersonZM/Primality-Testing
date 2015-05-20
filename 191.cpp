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
    int t;
    bool inter;
    Segmento s, s1, s2, s3, s4;
    Ponto p1, p2, p3, p4, p5, p6;
    int xleft, xright, ytop, ybottom;
    cin >> t;
    while(t--){
        cin >> p1.x >> p1.y >> p2.x >> p2.y >> xleft >> ytop >> xright >> ybottom;

        if(ytop < ybottom) swap(ytop,ybottom);
        if(xright < xleft) swap(xright, xleft);

        p3.x = xleft; p3.y = ybottom;
        p4.x = xleft; p4.y = ytop;
        p5.x = xright; p5.y = ytop;
        p6.x = xright; p6.y = ybottom;

        s.p = p1; s.q = p2;
        s1.p = p4; s1.q = p3;
        s2.p = p3; s2.q = p6;
        s3.p = p6; s3.q = p5;
        s4.p = p5; s4.q = p4;

        inter = false;
        if((p1.x >= xleft) && (p1.x <= xright) && (p1.y >= ybottom) && (p1.y <= ytop) ) inter = true;
        if((p2.x >= xleft) && (p2.x <= xright) && (p2.y >= ybottom) && (p2.y <= ytop) ) inter = true;
        if(Intercepta(s,s1)) inter = true;
        if(Intercepta(s,s2)) inter = true;
        if(Intercepta(s,s3)) inter = true;
        if(Intercepta(s,s4)) inter = true;
        if( (p1.x == p3.x && p1.y == p3.y) || (p1.x == p4.x && p1.y == p4.y) || (p1.x == p5.x && p1.y == p5.y) || (p1.x == p6.x && p1.y == p6.y)) inter = true;
        if( (p2.x == p3.x && p2.y == p3.y) || (p2.x == p4.x && p2.y == p4.y) || (p2.x == p5.x && p2.y == p5.y) || (p2.x == p6.x && p2.y == p6.y)) inter = true;

        if(inter) cout << "T" << endl;
        else cout << "F" << endl;
    }
   return 0;
}
