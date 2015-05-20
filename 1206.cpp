//Anderson Zudio de Moraes
//Julio Cesar
//Victor Cracel Messner
//AC 07/06/2014
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sstream>
using namespace std;
typedef struct {double x, y;}Ponto;
Ponto pol[2000000], pols[2000000]; //vetores globais

int SentidoPercurso(Ponto p1, Ponto p2, Ponto p3)
{  double a, b, c, d;
   a= p2.x - p1.x;  b= p3.y - p1.y;  c=p3.x - p1.x;  d=p2.y - p1.y;
   a = a*b-c*d;
   if (a > 0) return 1;
   else if (a < 0) return -1;
   else return 0;
}
void Quicksort(int e, int d)
// O ponto 0 fica fixo na ordenação, pois tem o menor x
{  int i,j,k;  Ponto t,r;  double d0, d1, d2;
   if (d > e)
   {  i=e;  j=d;  t=pol[(e+d)/2];
      d1= pol[0].x-t.x;  d2= pol[0].y-t.y;  d0=d1*d1+d2*d2;
      while(i<=j)
      {  while(1)
         {  k=SentidoPercurso(pol[0],pol[i],t);
            if (k==1) i++;
            else if (k==0)
            {  d1= pol[0].x-pol[i].x;  d2= pol[0].y-pol[i].y;  d1= d1*d1+d2*d2;
               if (d1 < d0) i++;
               else break;
            }
            else break;
         }
         while(1)
         {  k=SentidoPercurso(pol[0],t,pol[j]);
            if (k==1) j--;
            else if (k==0)
            {  d1= pol[0].x-pol[j].x;  d2= pol[0].y-pol[j].y;  d1= d1*d1+d2*d2;
               if (d0 < d1) j--;
               else break;
            }
            else break;
         }
         if (i<=j) {r=pol[i];  pol[i]=pol[j];  pol[j]=r;  i++; j--;}
      }
      Quicksort(e,j);  Quicksort(i,d);
   }
}
void OrdenaPontosy(int tp)
/*   Ordena pontos por ângulo, escolhendo como referência o ponto mais embaixo e mais
   à esquerda. */
{  int i,m;  Ponto p;
   m=0;
   for (i=1;i<tp;i++)
      if ((pol[i].y < pol[m].y)||(pol[i].y==pol[m].y)&&(pol[i].x<pol[m].x)) m=i;
   p=pol[0];  pol[0]=pol[m];  pol[m]=p;
   Quicksort(1,tp-1);
}
void OrdenaPontosx(int tp)
/* Ordena pontos por ângulo, escolhendo como referência o ponto mais à esquerda e mais
   embaixo. */
{  int i,m;  Ponto p;
   m=0;
   for (i=1;i<tp;i++)
      if ((pol[i].x < pol[m].x)||(pol[i].x==pol[m].x)&&(pol[i].y<pol[m].y)) m=i;
   p=pol[0];  pol[0]=pol[m];  pol[m]=p;
   Quicksort(1,tp-1);
}

int ConvexHull(int t)
{   int i,j,topo;
    OrdenaPontosy(t);
    /*Repete o ponto inicial no final*/
    pols[0]= pol[0]; pol[t]= pol[0];
    j= 1;
    /* Elimina pontos iniciais colineares */
    while (SentidoPercurso(pols[0],pol[j],pol[j+1])==0) j++;
    pols[1]= pol[j];  pols[2]= pol[j+1];
    topo= 2;
    for (i= j+2; i<= t;i++)
    {   while (SentidoPercurso(pols[topo-1], pols[topo], pol[i]) <= 0) topo--;
        topo++;  pols[topo]= pol[i];
    }
    return (topo);
}

void ExemploContornoConvexo(int np){
    int i, nv;
    nv = ConvexHull(np);
    for (i=0; i<=nv; i++) cout<< pols[i].x <<" " <<pols[i].y <<" ";
    cout<<endl<<endl;
}

int main()
{
    istringstream ss;
    string s;
    char lixo;
    int i, nv;
    double x, y;

    while(getline(cin, s)){
        ss.clear();
        ss.str(s);
        i = 0;
        while(ss >> lixo >> x >> lixo >> y >> lixo){ pol[i].x  = x; pol[i++].y = y;}
        nv = ConvexHull(i);
         for (i=0; i<nv; i++)  printf("(%g,%g) ",pols[i].x,pols[i].y);
		 printf("(%g,%g)\n",pols[i].x,pols[i].y);
    }


   return 0;
}
