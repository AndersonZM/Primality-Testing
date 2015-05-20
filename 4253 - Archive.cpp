//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//--WA
// -- Agradecimento especial para Thiago Gomes de Sá Farias, Giancarlo França e Raquel Marcolino --

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#define QZERO 0.000001
#define QZERO2 0.000001
//const double PI = acos(-1);

typedef struct Segmento{
    double y;
    double x1;
    double x2;
}segmento;

segmento S[5010];

using namespace std;


bool compare(segmento a, segmento b){
    return a.y < b.y;
}

int main(){
    int t, N;
    bool achei;
    double W, ang1, ang2, arq, d1, d2;
    pair <double, double> range, cone;

    scanf("%d",&t);
    while(t--){
        scanf("%lf %d", &W, &N);
        for(int i = 0; i < N; i++)
            scanf("%lf %lf %lf", &S[i].y, &S[i].x1, &S[i].x2);

        sort(S, S+N, compare);
        range = pair <double, double> (0,W);
        achei = false;

        //for(int i = 0; i < N; i++)
         //   cout << S[i].y << " " << S[i].x1 << " " << S[i].x2 << endl;


        while(arq = (range.second + range.first)/2, (range.second - range.first > QZERO) && !achei){
                d1 = S[0].x1 - arq;
                d2 = S[0].x2 - arq;

                //cout << d1 << " " << d2 << " " << d3 << " " << arq << endl;
               // cout << "Arq: " << arq << " " << S[0].x1 << " " << d2 << endl;
                ang1 = atan2(S[0].y,d2);
                ang2 = atan2(S[0].y,d1);
                cone = pair<double, double> (ang1,ang2);
                //cout << "Prim Segmento: " << S[0].y << " " << S[0].x1 << " " << S[0].x2 << endl;

                //cout << "Ang1 e Ang2: " <<(ang1*180)/PI << " " << (ang2*180)/PI << endl;

                for(int i = 1; i < N; i++){
                    d1 = S[i].x1 - arq;
                    d2 = S[i].x2 - arq;

                    ang1 = atan2(S[i].y,d2);
                    ang2 = atan2(S[i].y,d1);

                    //cout << "Cone1 e Cone2 Ang1 e Ang2: " << (cone.first*180)/PI << " " << (cone.second*180)/PI << " " << (ang1*180)/PI << " " << (ang2*180)/PI << endl;


                    if (ang1 - cone.second > QZERO){
                            range.first = arq;
                            break;
                    }
                    else if(cone.first - ang2 > QZERO){
                        range.second = arq;
                        break;
                    }


                    if(i == N-1){
                        achei = true;
                        //cout << "ultim Segmento " << S[i].y;
                        break;
                    }
                    if(ang1 - cone.first > QZERO2) cone.first = ang1;
                    if(cone.second - ang2 > QZERO2) cone.second = ang2;
                }
                //cout << (ang1*180)/PI << " " << (ang2*180)/PI << endl;
        }

        achei ? printf("YES\n") : printf("NO\n");
    }
    return 0;
}
