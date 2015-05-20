//Anderson Zudio
//Julio Cesar
//Victor Cracel messner
//AC em 07/06/2014
#include<iostream>
#include<algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
typedef struct {int x, y; char meudeus;}Ponto;
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
double AreaPoligono(int n)
{  long long c;  int i;
   c = 0;  for (i=0; i<n; i++) c= c + pol[i].x*pol[i+1].y-pol[i+1].x*pol[i].y;
   return c/2.0;
}
void Quicksort(int e, int d)
// O ponto 0 fica fixo na ordenação, pois tem o menor x
{  int i,j,k;  Ponto t,r;  long long d0, d1, d2;
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
void CaminhoFechado(int tp)
/* Ordena pontos e, ao final da ordenação, se houver pontos colineares com
   o ponto inicial, no final do vetor, a ordem dos pontos colineares é invertida. */
{  int i,m;  Ponto p;
   OrdenaPontosx(tp);
   i = tp-1;  while ((i>0)&&(SentidoPercurso(pol[0],pol[i-1],pol[i])==0)) i--;
   for (m=i; m<(i+tp)/2;m++) {p=pol[m]; pol[m]=pol[tp-1-m+i]; pol[tp-1-m+i]=p;}
   
   
   
}

int ConvexHull(int t)
{   int i,j,topo;
    OrdenaPontosx(t);
    /*Repete o ponto inicial no final*/
    pols[0]= pol[0]; pol[t]= pol[0];
    j= 1;
    pols[1]= pol[j];  pols[2]= pol[j+1];
    topo= 2;
    for (i= j+2; i<= t;i++)
    {   //while ((SentidoPercurso(pols[topo-1], pols[topo], pol[i]) <= 0) ) topo--;
        topo++;  pols[topo] = pol[i];
    }
    return (topo);
}

bool cmp(Ponto a, Ponto b){
		if(a.x == b.x)
			return a.y > b.y;
		else
			return a.x > b.x;
}

void ExemploContornoConvexo(){
    int i, np, nv;
    cin>> np;
    i = 0;
    while(i != np) {
        cin>> pol[i].x >>pol[i].y >> pol[i].meudeus;
        
        if(pol[i].meudeus == 'N')
        	np--;
        else
         	i++;
        
    }
    CaminhoFechado(np);    //nv = ConvexHull(np);
	//cout << nv << endl;
	cout << np << endl;
	//sort(pol,pol+np,cmp);
    for (i=0; i<np; i++) cout<< pol[i].x <<" " <<pol[i].y << endl;
}

int main()
{   int t;
	cin >> t;
	while(t--){
		ExemploContornoConvexo();
		
	}
   return 0;
}    
