//Anderson Zudio
//Julio Cesar
//Victor Cracel messner
//AC 24/06/2014
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define QZERO 0.000000001
#define PI acos(-1)
typedef struct {double x, y;}Ponto;
typedef struct segmento{Ponto p, q;} Segmento;
Ponto pol[100000], pols[100000];

using namespace std;

/* Dado um conjunto de pontos que formam um polígono qualquer, o seu contorno convexo
tem pelo menos um dos lados coincidindo com um lado do seu menor retangulo que cobre todo o polígono.

A partir dessa definição, o algoritmo a seguir é uma implementação que descobre quais dos
lados do contorno do convexo faz o menor retangulo que cobre o polígono calculando somente a área dos retangulos*/


//Usado na interceção de segmentos
int SentidoPercurso(Ponto p1, Ponto p2, Ponto p3)
{  double a, b, c, d;
   a= p2.x - p1.x;  b= p3.y - p1.y;  c=p3.x - p1.x;  d=p2.y - p1.y;
   a = a*b-c*d;
   if (a > 0) return 1;
   else if (a < 0) return -1;
   else return 0;
}

//Verifica a interceção de dois segmentos
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

//Usado no convex hull
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
//Usado pelo QuickSort
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

//Retorna o contorno do convexo de um poligono que está armazenado em pol
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




int main(){
    int n;
    cout.setf(ios::fixed);
    cout.precision(4);
    while(cin >> n, n){
            for(int i = 0; i < n; i++)
                cin >> pol[i].x >> pol[i].y;

            //Pegadinha do UVA, se n < 3 então não existe polígono. Então a área do retangulo é zero
            if(n < 3) {
                cout << "0.0000" << endl;
                continue;
            }

            //Contorno convexo
            auto nv = ConvexHull(n);

            //Possibilidade de não existir o polígono, o problema não garante nada sobre pontos colineares.
            if(nv < 3) {
                cout << "0.0000" << endl;
                continue;
            }

            //A variavel mArea guarda a a área minima do retangulo procurado.
            double mArea = -1;
            //O array pols fornece o contorno do convexo do polígono de entrada
            for(int i = 1; i <= nv; i++){
                Ponto vetor0, vetor1;
                vetor0.x = pols[i].x - pols[i-1].x; vetor0.y = pols[i].y - pols[i-1].y;
                double t = sqrt(vetor0.x*vetor0.x + vetor0.y*vetor0.y);
                vetor0.x /= t; vetor0.y /= t;
                vetor1.x = -vetor0.y; vetor1.y = vetor0.x;
                double s0, t0, s1, t1;
                s0 = t0 = s1 = t1 = 0;
                //Descobre os lados do retangulo alinhado com o segmento pols[i] e pols[i-1]
                for(int j = 1; j < nv; j++){
                    Ponto D;
                    D.x = pols[j].x - pols[0].x; D.y = pols[j].y - pols[0].y;
                    double test = vetor0.x*D.x + vetor0.y*D.y;
                    if(test < s0) s0 = test;
                    else if(test > s1) s1 = test;
                    test = vetor1.x*D.x + D.y*vetor1.y;
                    if(test < t0) t0 = test;
                    else if(test > t1) t1 = test;
                }
                //Calcula área do retangulo e guarda a menor
                double Area = (s1-s0)*(t1-t0);
                if(Area < mArea) mArea = Area;
                else if(mArea == -1) mArea = Area;
            }
            cout << mArea << endl;
    }
    return 0;
}
