//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//AC em 04/05/2014 

#include <iostream>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <vector>
#define MAX 2000
 
using namespace std;
 
int main() {
    int T;
    vector <int> S_1, S_2;
    int n, i, j, seq[MAX], s_max = 0, g_max = 0;
    cin >> T;
    while(T--) {
        memset(seq, 0, sizeof(seq));
        cin >> n;
        for (i = 0; i < n; i++)
            S_1.push_back((cin >>  j,j));
        cin >> n;
        for (i = 0; i < n; i++)
            S_2.push_back((cin >>  j,j));
 
        for (i = 0; i < S_1.size(); i++) {
            s_max = 0;
            for (j = 0; j < S_2.size(); j++) {
                if (S_1.at(i) > S_2.at(j) && s_max < seq[j]) s_max = seq[j];
                if (S_1.at(i) == S_2.at(j)) seq[j] = (1 + s_max);
                g_max = max(g_max,max(s_max, seq[j]));
            }
        }
 
        cout << g_max << endl;
        g_max = 0;
        S_1.clear();
        S_2.clear();
    }
    return 0;
}