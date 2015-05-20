//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
//AC em 10/05/2014

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#define MAX 11000
#define QUERY 1100
#define MODULO(a) ((a) >= 0) ? (a) : (-a)
using namespace std;
 
struct dicionario {
    char str[30];
    int t;
    int next;
};
 
int main() {
    dicionario D[MAX], W;
    int n, q, i=1, j, pos_i[30], next[30];
    int a, b, c, n_e, r_e, erros = 0;
    char str[30];
 
    cin >> n;
    for (j = 0; j < sizeof(pos_i); j++) pos_i[j] = n+1;
    memset(next,0,sizeof(next));
    while(n >= i) {
            cin >> D[i].str;
            D[i].t = strlen(D[i].str);
            D[i].next = (n+1);
 
            if (!next[D[i].t]) {
               pos_i[D[i].t] = i;
               next[D[i].t] = i;
            }
            else {
               D[next[D[i].t]].next = i;
               next[D[i].t] = i;
            }
    i++;
    }
    cin >> q;
    while(q--) {
          cin >> str;
          r_e = (n+1);
          a = pos_i[strlen(str)-1];
          b = pos_i[strlen(str)];
          c = pos_i[(strlen(str)+1)];
          while (n_e = min(c,min(a,b)), n_e <= n) {
            W = D[n_e];
            erros = 0;
            if (strcmp (W.str,str) != 0){
                if (W.t != strlen(str)) {
                    if (MODULO((W.t-strlen(str))) == 1) {
                        i = j = 0;
                        while(i <= W.t && j <= strlen(str) && erros <= 1)
                            if (W.str[j] != str[i]) {
                                erros++;
                                if (W.t > strlen(str) ) j++;
                                else i++;
                            }
                            else { j++; i++; }
                    }
                    (n_e == a) ? a = D[n_e].next : c = D[n_e].next;
                }
                else {
                    for (j=0; j<W.t; j++)
                        if (W.str[j] != str[j]) {
                            erros++;
                            if (erros>1) break;
                            else
                                if (W.str[(j+1)] == str[j] && W.str[j] == str[(j+1)]) j++;
 
                        }
                    b = D[n_e].next;
                }
            }
            else {
                r_e = n_e;
                break;
            }
 
            if (erros == 1) {
                r_e = min(r_e,n_e);
                a = n+1;
                c = n+1;
            }
 
          }
 
          if (r_e <= n) {
            if (strcmp (D[r_e].str,str) != 0)
                cout << str << " is a misspelling of " << D[r_e].str << endl;
            else
                cout << str << " is correct" << endl;
          }
          else
            cout << str << " is unknown" << endl;
    }
 
    return 0;
}