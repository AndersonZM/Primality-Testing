//Anderson Zudio
//AC em 17/06/2014
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
typedef struct {int x, y;}Ponto;
typedef struct {Ponto p, q;}Segmento;
Ponto pol[505], pols[505], pols2[505]; //vetores globais
int M, C, Mv, Cv;
bool eC, eM;

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

int ConvexHull(int t, bool quem)
{   int i,j,topo;
    OrdenaPontosy(t);
    /*Repete o ponto inicial no final*/
    pols[0] = pol[0]; pol[t] = pol[0];
    j= 1;
    /* Elimina pontos iniciais colineares */
    while (SentidoPercurso(pols[0],pol[j],pol[j+1]) == 0 && j <= t) j++;
    pols[1] = pol[j];  pols[2]= pol[j+1];
    topo = 2;
    if(quem){
        if(j > t){
            for(i = 0; i <= t; i++)
                pols[i] = pol[i];
            eC = false;
            return t;
        }
    }
    else{
        if(j > t){
            for(i = 0; i <= t; i++)
                pols[i] = pol[i];
            eM = false;
            return t;
        }
    }
    for (i= j+2; i<= t;i++)
    {   while (SentidoPercurso(pols[topo-1], pols[topo], pol[i]) <= 0) topo--;
        topo++;  pols[topo] = pol[i];
    }
    return (topo);
}

void ExemploSentidoPercurso(){
    Ponto p1,p2,p3;  int j;
    cout<<"Exemplo de uso de SentidoPercurso"<<endl;
    cout<< "Entre com as coordenadas (x,y) de p1, p2 e p3: ";
    cin>>p1.x >>p1.y >>p2.x >>p2.y >>p3.x >>p3.y;
    j= SentidoPercurso(p1,p2,p3);
    if (j==1)       cout<<"Antihorario"<<endl;
    else if (j==-1) cout<<"Horario" <<endl;
    else            cout<<"Pontos colineares"<<endl;
    cout<<endl;
}

void ExemploIntersecaoSegmentos(){
    Segmento s1,s2;
    cout<<"Exemplo de uso de Intercepta" <<endl;
    cout<<"Entre com o segmento 1 (p1 e p2): ";
    cin>> s1.p.x >> s1.p.y >> s1.q.x >> s1.q.y;
    cout<< "Entre com o segmento 2 (p3 e p4): ";
    cin>> s2.p.x >> s2.p.y >> s2.q.x >> s2.q.y;
    if (Intercepta(s1,s2)==1) cout<< "Interceptam" <<endl;
    else cout<< "Nao interceptam" <<endl;
    cout<<endl;
}


bool PontoInterior(Ponto *mpol ,Ponto q, int n){
    Segmento s, s1;
    Ponto qf;
    int cont;
    mpol[0] = mpol[n];
    for(int i = 1; i <= n; i++){
        s.p = mpol[i-1]; s.q = mpol[i]; s1.p = q; s1.q = q;
        if(((mpol[i].x == q.x) && (mpol[i].y == q.y)) || Intercepta(s,s1)) return true;
    }
    qf.x = 99999; qf.y = q.y; cont = 0;

    for(int i = 1; i <= n; i++){
        s.p = q; s.q = qf; s1.p = mpol[i-1]; s1.q = mpol[i];
        if((((mpol[i].y > q.y) && (mpol[i-1].y <= q.y)) || ((mpol[i-1].y > q.y) && (mpol[i].y <= q.y))) && (Intercepta(s,s1)))
            cont++;
    }
    return((cont % 2) != 0);


}

bool resolve(){
    Segmento s1, s2;
    //ESPECIAL
    if(Mv == 1 && Cv == 1){
        if(pols2[0].x == pols[0].x && pols2[0].y == pols[0].y) return false;
        else return true;
    }
    else if(Mv == 1 && Cv == 2){
        s1.p = pols[0]; s1.q = pols[1]; s2.p = s2.q = pols2[0];
        if(Intercepta(s1, s2)) return false;
        else return true;
    }
    else if(Mv == 2 && Cv == 1){
        s1.p = pols2[0]; s1.q = pols2[1]; s2.p = s2.q = pols[0];
        if(Intercepta(s1, s2)) return false;
        else return true;
    }
    else if(Mv == 2 && Cv == 2){
        s1.p = pols[0]; s1.q = pols[1]; s2.p = pols2[0]; s2.q = pols2[1];
        if(Intercepta(s1, s2)) return false;
        else return true;
    }
    //FIM CASOS ESPECIAIS
    if(eM)
        for(auto i = 0; i < Cv; i++)
            if(PontoInterior(pols2, pols[i], Mv)) return false;
    if(eC)
        for(auto i = 0; i < Mv; i++)
            if(PontoInterior(pols, pols2[i], Cv)) return false;

    for(auto i = 0; i < Mv; i++){
        s1.p = pols2[i]; s1.q = pols2[i+1];
        for(auto j = 0; j < Cv; j++){
            s2.p = pols[j]; s2.q = pols[j+1];
            if(Intercepta(s1, s2)) return false;
        }
    }
    return true;

}

int main()
{
   while(cin >> M >> C, M){
        eM = eC = true;
        for(auto i = 0; i < M; i++)
            cin >> pol[i].x >> pol[i].y;

        if(M > 2)
            Mv = ConvexHull(M, false);
        else{
            Mv = M;
            eM = false;
            for(auto i = 0; i < Mv; i++)
                pols[i] = pol[i];
            pols[Mv] = pols[0];
        }

        for(auto i = 0; i <= Mv; i++){
            pols2[i] = pols[i];
        }

        for(auto i = 0; i < C; i++)
            cin >> pol[i].x >> pol[i].y;

        if(C > 2)
            Cv = ConvexHull(C, true);
        else{
            Cv = C;
            eC = false;
            for(auto i = 0; i < Cv; i++)
                pols[i] = pol[i];
            pols[Cv] = pols[0];
        }

        resolve() ? cout << "Yes" : cout << "No";
        cout << endl;

        /*for(auto i = 0; i < Mv; i++)
            cout << pols2[i].x << " " << pols2[i].y << " ";
        cout << endl;

        for(auto i = 0; i < Cv; i++)
            cout << pols[i].x << " " << pols[i].y << " ";
        cout << endl;*/


   }
   return 0;
}
