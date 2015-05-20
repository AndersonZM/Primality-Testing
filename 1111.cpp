//Anderson Zudio
//AC em 22/07/2014
#include<iostream>
#include <complex>
#include <math.h>
#define inf 1001000
#define QZERO 0.000000001
#define PI acos(-1)

using namespace std;
typedef struct {double x, y;}Ponto;
typedef struct {Ponto p, q;}Segmento;
typedef struct { double a,b,c; }reta;
Ponto pol[100000], pols[100000]; //vetores globais

int SentidoPercurso(Ponto p1, Ponto p2, Ponto p3)
{  long long a, b, c, d;
   a= p2.x - p1.x;  b= p3.y - p1.y;  c=p3.x - p1.x;  d=p2.y - p1.y;
   a = a*b-c*d;
   if (a > 0) return 1;
   else if (a < 0) return -1;
   else return 0;
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

double RotatingCalipers(int n){
    int a, b;
    double angulo_r = 0, mdia, miny, maxy, anguloa, angulob;
    complex<double> caliper_a, caliper_b, edge_a, edge_b, aux_a, aux_b;
    caliper_a = complex<double>(1,0); caliper_b = complex<double>(-1,0);

    miny = maxy = inf; maxy = -maxy;
    for(int i = 0; i < n; i++){
        if(maxy < pols[i].y){
            maxy = pols[i].y;
            b = i;
        }
        if(miny > pols[i].y){
            miny = pols[i].y;
            a = i;
        }
    }
    mdia = pols[b].y - pols[a].y;
    while(angulo_r < PI){
        edge_a  = complex<double>(pols[a + 1].x - pols[a].x, pols[a + 1].y - pols[a].y);
        edge_b  = complex<double>(pols[b + 1].x - pols[b].x, pols[b + 1].y - pols[b].y);
        aux_a   = complex<double>(pols[a].x, pols[a].y);
        aux_b   = complex<double>(pols[b].x, pols[b].y);
        anguloa = acos( ((conj(edge_a)*caliper_a).real()) / abs(caliper_a) / abs(edge_a));
        angulob = acos( ((conj(edge_b)*caliper_b).real()) / abs(caliper_b) / abs(edge_b));
        caliper_a = caliper_a * polar((double) 1, min(anguloa, angulob));
        caliper_b = caliper_b * polar((double) 1, min(anguloa, angulob));

        if(anguloa < angulob){
            mdia = min(mdia,(abs(((conj((aux_a + caliper_a)-(aux_a))*((aux_b)-(aux_a))).imag())/abs((aux_a + caliper_a)-(aux_a)))));
            a = (a + 1) % n;
        }
        else{
            mdia = min(mdia,(abs(((conj((aux_b + caliper_b)-(aux_b))*((aux_a)-(aux_b))).imag())/abs((aux_b + caliper_b)-(aux_b)))));
            b = (b + 1) % n;
        }
        angulo_r += min(anguloa, angulob);
    }
    return mdia;
}

int main(){
    double resp;
    int n, t, nv;
    cout.setf(ios::fixed);
    cout.precision(2);
    t = 0;
    while(cin >> n, n){
        for(int i = 0; i < n; i++)
            cin >> pol[i].x >>pol[i].y;


        nv = ConvexHull(n);

        resp = RotatingCalipers(nv);
        cout << "Case " << ++t << ": " << resp << endl;

    }
   return 0;
}
