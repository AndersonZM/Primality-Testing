//Anderson Zudio de Moraes
//Julio Cesar
//Victor Cracel Messner
//AC 31/05/2014
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <deque>

using namespace std;
typedef struct {double x, y;}Ponto;
typedef struct {Ponto p, q;}Segmento;

bool inter[100010];
Segmento s[100010];
int i;
deque <int> fila;

int SentidoPercurso(Ponto p1, Ponto p2, Ponto p3)
{  double a, b, c, d;
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
    int t, m, aux, n;
    while(cin >> m, m){

        for(i = 1; i <= m; i++){
            cin >> s[i].p.x >> s[i].p.y >> s[i].q.x >> s[i].q.y;
            fila.push_back(i);

            n = fila.size();
            for(int j = n-2; j >= 0; j--){
                    if(Intercepta(s[i],s[fila[j]])){
                        fila.erase(fila.begin()+j);
                    }
            }


            /*for(int x = 0; x < fila.size(); x++)
                cout << fila[x] << " ";
            cout << endl;*/

        }
    cout << "Top sticks: ";

    while(!fila.empty()){
        aux = fila.front();
        fila.pop_front();
        cout << aux;
        if(fila.empty()) cout << "." << endl;
        else cout << ", ";
    }

}
   return 0;
}
