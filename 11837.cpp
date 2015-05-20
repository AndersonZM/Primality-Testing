//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar 
//AC em 22/05/2014

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#define MODMUSICAL(a,b) ((a) >= (b)) ? (1+(a-b)) : (13-(b-a))
using namespace std;
 
char retNota(char *nota, int len) {
    int inc = 0;
    switch (nota[0]-64){
        case 7: inc++;
        case 6: case 5: inc++;
        case 4: inc++;
        case 3: case 2: inc++;
        default: break;
    }
    return (nota[0] + inc + ((len == 2) ? ((nota[1] == '#') ? 1:(-1)) : 0));
}
 
int main() {
    char notas[100010], difs[10010];
    int M,T;
    char str[3], nt, ant, inicio;
    int a,p, i,j;
    bool plagio;
 
  while(scanf("%i %i", &M, &T), M && T) {
        plagio = false;
 
        scanf("%s", str);
        ant = ( nt = retNota(str, strlen(str)), (nt) >= 65 ? nt : 76);
        for (i = 1; i < M ; i++) {
            scanf("%s", str);
            nt = retNota(str, strlen(str));
            if (nt == 64) nt = 76;
            a = MODMUSICAL(nt,ant);
            notas[(i-1)] = 64 + a;
            ant = nt;
        }
 
        scanf("%s", str);
        ant = ( nt = retNota(str, strlen(str)), (nt) >= 65 ? nt : 76);
        for (i = 1; i < T ; i++) {
            scanf("%s", str);
            nt = retNota(str, strlen(str));
            if (nt == 64) nt = 76;
            a = MODMUSICAL(nt,ant);
            difs[(i-1)] = 64 + a;
            ant = nt;
        }
 
        if (T > 1) {
            string Otrecho(notas); string Oplag(difs); string plag;
            inicio = Oplag[0];
            for (i = 0; i < (M-1) ; i++)
                if (Otrecho[i] == inicio) {
                    if (M < (i+(T-1)) ) break;
                    plag = Otrecho.substr (i,(T-1));
                    if ( Oplag == plag ) {
                        plagio = true;
                        break;
                    }
                }
        }
 
        if (plagio) printf("S\n");
        else printf("N\n");
    memset(difs,0,T*sizeof(char));
    memset(notas,0,M*sizeof(char));
    }
    return 0;
}
