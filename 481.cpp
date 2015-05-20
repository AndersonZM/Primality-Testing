//AC em 26/04/2014
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>

using namespace std;

int T[10000001], S[10000001], O[10000001], n, m, kscm; //TM[10001][10001];

int SCM() {
    int j, k = 1;
	T[1] = S[1];   O[1] = 1;
	for(int i = 2; i <= n; i++ ) {
        if( S[i] > T[k] ) {
            T[++k] = S[i];    O[i] = k;
        } else {
            j = lower_bound( T, T + k, S[i] ) - T;
            if( S[i] > T[j] )  j++;
            T[j] = S[i];  O[i] = j;
        }
    }
    return k;
}

void ImprimeSCM(int i,int j,int m ) {
     if( j > 0 ) {
         if( ( O[i] == j ) && ( S[i] < m ) ) {
             ImprimeSCM( i-1, j-1, S[i] );
             cout << S[i] << endl;
         }
         else ImprimeSCM( i-1, j, m );
     }
}

int main()
{
    n = 0;
    while(cin >> S[++n]);
    kscm = SCM();
    cout << kscm << endl << '-' << endl;
    ImprimeSCM( n, kscm, 2147483647);
    return 0;
}
