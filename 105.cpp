//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//AC em 08/06/2014

#include <string.h>
#include <algorithm>
#include <cstdio>

using namespace std;

int main() {
    int pred[10001];
    int e, h, d, i, pos;
    bool f=true;
    memset(pred,0,10000*sizeof(int));

	while (scanf("%d %d %d", &e, &h, &d) != EOF)
		for (i = e; i < d; i++)
            if (h > pred[i]) pred[i] = h;


    i = 0;
    for (pos = 1; pos < 10000; pos++)
        if (pred[pos] != i) {
           (!f) ? printf(" ") : (f=false);
           printf("%d %d", pos, pred[pos]);
           i = pred[pos];
        }
    printf("\n");

	return 0;
}