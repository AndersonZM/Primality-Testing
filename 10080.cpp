//Anderson Zudio
//AC 28/07/2014
#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#define INF 0x3f3f3f3f
#define QZERO 0.00000001
#define NVM 1006
#include <stdlib.h>
using namespace std;

struct nof{int v; int p; int ma;};

int C[NVM][NVM], CC[NVM][NVM], F[NVM][NVM], FC[NVM][NVM], n, m;
bool V[NVM];  nof Fila[NVM];

void Imprime(){
        cout<<"Fluxo:"<<endl;
        for(int i=1;i<=n;i++) {
            for(int j=1;j<=n;j++) {
                if (F[i][j] > 0) cout << i << " " << j << "  "<<F[i][j]<<endl;
            }
        }
}
void ImprimeC(){
        cout<<"Rede:"<<endl;
        for(int i=1;i<=n;i++) {
            for(int j=1;j<=n;j++) {
                if (CC[i][j] > 0) cout << i << " " << j << "  "<<CC[i][j]<<endl;
            }
        }
}
int MFBFS(){
    int i, j, f, r, p, s, fm;
    memset(V, 0, sizeof(V));
    f=r=1;   Fila[1].v=1; Fila[1].p=0;  Fila[1].ma=INF; V[1]=true;
    while(f<=r){
        s=Fila[f].v;
        for(int i=1; i<=n; i++) if (!V[i]){
            if (C[s][i] > 0){
                V[i]=true;
                Fila[++r].v=i;  Fila[r].p=f;  Fila[r].ma=min(Fila[f].ma,C[s][i]);
                if (i==n) {f=r; break;}
            }
        }
        f++;
    }
    if (Fila[r].v==n){
        j=r;  i=Fila[j].v;  p=Fila[j].p;  s=Fila[p].v;  fm=Fila[r].ma;
        while (i != 1){
            C[s][i] -= fm;  C[i][s] += fm;
            F[s][i] += fm;  F[i][s] -= fm;
            j=p;  i=Fila[j].v;  p=Fila[j].p;
            if (i != 1) s= Fila[p].v;
        }
        return fm;
    }
    return 0;
}

double dist(pair <double, double> p, pair <double, double> q){
    return ((p.first-q.first)*(p.first-q.first)) + ((p.second-q.second)*(p.second-q.second));
}

int main() {
    int fmax;
    int mt,u,nv,nva,cp,ck;
    int N, M, S, V;
    double x, y;
    vector <pair <double, double> > gophers, holes;
    n=3;
    while(cin >> N >> M >> S >> V){
        gophers.clear();
        holes.clear();
        for(int i = 0; i < N; i++){
            cin >> x >> y;
            gophers.push_back({x, y});
        }
        for(int i = 0; i < M; i++){
            cin >> x >> y;
            holes.push_back({x, y });
        }

       n = N + M + 2;
       m = 0;
       for(int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++){
             if(i == 1 && j > 1 && j <= N+1){
                 CC[i][j] = 1;
                 m++;
             }
             else if(i == n && j > N+1 && j <= N+M+1){
                 CC[j][i] = 1;
                 m++;
             }
             else{
                CC[i][j] = 0;
             }
             FC[i][j]=0;
        }
       double aux = V*S*S*V;
       for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
                if(dist(gophers[i], holes[j]) <= aux){
                    m++;
                    CC[i+2][N+j+2] = 1;
                }
            }
       }
       nv=20;
       nva=nv;
       while(nva--){
           for (int i=1; i<=n; i++) for (int j=1; j<=n; j++){
               C[i][j]=CC[i][j];   F[i][j]=FC[i][j];
           }
           fmax = 0;
           while(int fm = MFBFS()) {
               fmax += fm;
           }
       }
       cout << N - fmax << endl;
    }
}
