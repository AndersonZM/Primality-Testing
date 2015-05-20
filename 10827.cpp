//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//AC em 02/05/2014

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;
int M[200][200], S[1000], n, n2;

void SCSM(){
    int sm, cm, fm, ss, cs, fs;
	sm = cm = fm = ss = cs = fs = 0;

	for( int i = 1; i <= n; i++ ) {
		if( ( ss + S[i] ) >= 0 ) {
			ss=ss + S[i];
			fs=i;
			if( cs == 0 ) cs=i;
			if( ss > sm ) {
                sm=ss;   cm=cs;   fm=fs;
            }
        }
		else {
			ss = cs = fs = 0;
        }
	}
	cout << endl << endl << "Soma Maxima = " << sm;
}

int main()
{
    int t, soma, asoma, asoma2;
    cin >> t;
    while(t--){
        cin >> n;
        for( int i = 1; i <= n; i++ )
            for(int j = 1; j <= n; j++){
                cin >> M[i][j];
                M[i + n][j] = M[i][j];
                M[i][j + n] = M[i][j];
                M[i + n][j + n] = M[i][j];
            }

        n2= 2*n;

        for(int i = 1; i <= n; i++)
            for(int j = 2; j <= n2; j++){
                M[i][j] += M[i][j - 1];
                M[i + n][j] = M[i][j];
            }

        soma = -10000000;
        for(int i = 1; i <= n2; i++){

            if(i > 1)
                for(int x = 1; x <= n2; x++)
                 for(int y = i; y <= n2; y++)
                    M[x][y] -= M[x][i-1];

            for(int k = i; k <= n2; k++){
                asoma = 0;
                for(int j = 1; j <= n2; j++){
                    asoma += M[j][k];
                    asoma2 = asoma;
                    for(int l = 1; l <= j; l++){
                        if((j - l + 1) <= n && (k - i + 1 <= n)){
                            soma = max(soma,asoma2);
                        }
                        asoma2 -= M[l][k];
                    }
                }
            }
        }

        cout << soma << endl;
        /*for(int i = 1; i <= n2; i++){
            for(int j = 1; j <= n2; j++)
                cout << M[i][j] << " ";
            cout << endl;
        }*/
    }



    return 0;
}
