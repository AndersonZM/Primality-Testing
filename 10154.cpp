//AC 01/05/2014
//Julio Cesar
//Anderson Zudio
//Victor Cracel

#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#define MAX 6000
const int TETO = 0x7ffffff;

using namespace std;

struct Turtle { int p,c;
    Turtle(int p, int c) : p(p), c(c) {}
};

bool comparar(const Turtle& a, const Turtle& b) {
    return (a.c != b.c) ? (a.c < b.c) : (a.p > b.p);
}

int sol[MAX];
int main() {
    vector <Turtle> V;
    int p, c, p_t, tamanho=0;
    while(cin >> p >> c)
        if ((c-p) >= 0)
            V.push_back(Turtle(p, c));

    sort(V.begin(), V.end(), comparar);
    if (V.size() > 2 && ((V[0].p + V.back().p) <= V.back().c)) {
        for (int i = 1; i <= V.size(); i++) sol[i] = TETO;
            for(int i = 0; i <= V.size(); i++)
                for(int j = tamanho; j >= 0; j--) {
                    if ( (sol[j] + V[i].p <= V[i].c) && (sol[j] + V[i].p < sol[j + 1]) ) {
                        sol[j + 1] = sol[j] + V[i].p;
                        tamanho = max(j + 1,tamanho);
                    }
                }
    } else tamanho = (V.size() > 1) ? (1 + (int) ((V[0].p + V.back().p) <= V.back().c) ) : V.size();
    cout << tamanho << endl;
    return 0;
}
