//Anderson Zudio
//Julio Cesar
//Victor Cracel
//AC em 31/05/2014

#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ABS(a) (a)>0 ? (a):(-a)
using namespace std;
 
struct ponto {
    long long x,y;
    ponto(long long x,long long y) : x(x), y(y) {}
    ponto() {}
};
 
typedef struct {ponto p, q;} Segmento;
ponto pol[100000], pols[100000]; //vetores globais
 
int SentidoPercurso(ponto p1, ponto p2, ponto p3)
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
{  double c;  int i,j;
   c = 0;  
   for (i=0; i<n; i++) {
   j = (i+1)%n;
   c += (pol[i].x+pol[j].x)*(pol[i].y-pol[j].y);
   }
   return c/2;
}
 
double AreaConvexa(int n)
{  long long c;  int i;
   c = 0;  for (i=0; i<n; i++) c= c + pols[i].x*pols[i+1].y-pols[i+1].x*pols[i].y;
   return c/2.0;
}
void Quicksort(int e, int d)
// O ponto 0 fica fixo na ordenação, pois tem o menor x
{  int i,j,k;  ponto t,r;  long long d0, d1, d2;
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
 
{  int i,m;  ponto p;
   m=0;
   for (i=1;i<tp;i++)
      if ((pol[i].y < pol[m].y)||(pol[i].y==pol[m].y)&&(pol[i].x<pol[m].x)) m=i;
   p=pol[0];  pol[0]=pol[m];  pol[m]=p;
   Quicksort(1,tp-1);
}    
void OrdenaPontosx(int tp)
 
{  int i,m;  ponto p;
   m=0;
   for (i=1;i<tp;i++)
      if ((pol[i].x < pol[m].x)||(pol[i].x==pol[m].x)&&(pol[i].y<pol[m].y)) m=i;
   p=pol[0];  pol[0]=pol[m];  pol[m]=p;
   Quicksort(1,tp-1);
}
void CaminhoFechado(int tp)
 
{  int i,m;  ponto p;
   OrdenaPontosx(tp);
   i = tp-1;  while ((i>0)&&(SentidoPercurso(pol[0],pol[i-1],pol[i])==0)) i--;
   for (m=i; m<(i+tp)/2;m++) {p=pol[m]; pol[m]=pol[tp-1-m+i]; pol[tp-1-m+i]=p;}
}
 
int ConvexHull(int t)
{   int i,j,topo;
    OrdenaPontosy(t);
   
    pols[0]= pol[0]; pol[t]= pol[0];
    j= 1;
 
    while (SentidoPercurso(pols[0],pol[j],pol[j+1])==0) j++;
    pols[1]= pol[j];  pols[2]= pol[j+1];
    topo= 2;
    for (i= j+2; i<= t;i++)
    {   while (SentidoPercurso(pols[topo-1], pols[topo], pol[i]) <= 0) topo--;
        topo++;  pols[topo]= pol[i];
    }
    return (topo);
}
 
int main() {  
    int N, i, a,b, I=0;
    double areaC, areaT, pct;
    cout.setf(ios_base::fixed);
    cout.precision(2);
    while(cin >> N, N){
        I++;
        for (i=0; i< N; i++)
           pol[i] = (cin >> a >> b, ponto(a,b));
       
        areaT = ABS(AreaPoligono(N));
        a = ConvexHull(N);
        areaC = ABS(AreaConvexa(a));
        pct = 100*(1-(areaT/areaC));
        cout << "Tile #" << I << endl << "Wasted Space = " << pct << " %" << endl << endl;
    }  
   return 0;
}