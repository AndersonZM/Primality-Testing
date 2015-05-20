//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
//Juan Lopes
// AC em 13/05/2014

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#define MAX 1500
using namespace std;

int RA[MAX], tempRA[MAX];
int SA[MAX], tempSA[MAX];
int C[MAX];
int Phi[MAX], PLCP[MAX], LCP[MAX], F[MAX];
char P[MAX], T[MAX];
int N, M, cont, R[256];

void BM()
{  int i,j,k;
    cont = 0;
    //Preprocessamento
    for(i = 0; i < 256; i++)  R[i] = M;
    for(i = 0; i < M-1; i++)  R[P[i]] = M-i-1;
    //Busca
    k = M-1;
    while(k < N) {
       i = k; j = M-1;
       while(T[i] == P[j] && j >= 0) {  i--;  j--;
       }
       if(j < 0) cont++;
       k+= R[T[k]];
    }
}

void suffix_sort(int n, int k) {
    memset(C, 0, sizeof C);

    for (int i = 0; i < n; i++)
        C[i + k < n ? RA[i + k] : 0]++;

    int sum = 0;
    for (int i = 0; i < max(256, n); i++) {
        int t = C[i];
        C[i] = sum;
        sum += t;
    }

    for (int i = 0; i < n; i++)
        tempSA[C[SA[i] + k < n ? RA[SA[i] + k] : 0]++] = SA[i];

    for (int i = 0; i < n; i++) SA[i] = tempSA[i];
}

void suffix_array(string &s) {
    int n = s.size(), r;

    for (int i = 0; i < n; i++){
        RA[i] = s[i] - 1;
        SA[i] = i;
    }

    for (int k = 1; k < n; k *= 2) {
        suffix_sort(n, k);
        suffix_sort(n, 0);

        r = tempRA[SA[0]] = 0;
        for (int i = 1; i < n; i++) {
            bool equal = true;
            equal &= RA[SA[i]] == RA[SA[i-1]];
            equal &= RA[SA[i]+k] == RA[SA[i-1]+k];

            tempRA[SA[i]] = equal ? r : ++r;
        }
        for (int i = 0; i < n; i++) RA[i] = tempRA[i];
    }
}

void lcp(string &s) {
    int n = s.size();

    Phi[SA[0]] = -1;
    for (int i = 1; i < n; i++)
        Phi[SA[i]] = SA[i-1];

    int L = 0;
    for (int i = 0; i < n; i++) {
        if (Phi[i] == -1) {
            PLCP[i] = 0;
            continue;
        }
        while (s[i + L] == s[Phi[i] + L])
            L++;

        PLCP[i] = L;
        L = max(L-1, 0);
    }

    for (int i = 1; i < n; i++)
        LCP[i] = PLCP[SA[i]];
}

int main()
 {
    int t, mlcp, ilcp, i, j, k;
    string s;
    cin >> t;
    while(t--){
        cin >> s;
        s += '.';
        suffix_array(s);
        lcp(s);

        mlcp = ilcp = 0;
        for (i = 0; i< s.length(); i++){
            if(LCP[i] > mlcp){
                ilcp = i;
                mlcp = LCP[i];
            }
        }

        if(!mlcp){
            cout << "No repetitions found!" << endl;
        }
        else{
            memset(P, '\0', MAX);
            memset(T, '\0', MAX);
            for(i = 0, j = SA[ilcp], k = 0; k < mlcp; j++, i++, k++){
                P[i] = s[j];
            }
            for(i = 0; i < s.length(); i++)
                T[i] = s[i];

            N = i-1;
            M = strlen(P);
            BM();
            cout << P << " " << cont << endl;
        }


        /*for(int i=1; i< s.length(); i++)
            cout << i;
        cout <<endl << endl;
        for (int i = 1; i< s.length(); i++)
            cout << SA[i]+1<<" ";
        cout <<endl;
        for (int i = 1; i< s.length(); i++)
            cout << LCP[i]<<" ";
        cout <<endl;*/

    }
    return 0;
}
