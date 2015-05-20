//AC em 26/04/2014
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>

using namespace std;
int n, m, TM[1010][1010];
char A[1010], B[10010];

void SMC(){
	n = strlen(A);
	m = strlen(B);
    for( int i=1; i<=n; i++ ) TM[i][0]=0;
    for( int j=1; j<=m; j++ ) TM[0][j]=0;
	for( int i=1; i<=n; i++ ) {
		for( int j=1; j<=m; j++ ) {
			if( A[i-1] == B[j-1] )
				TM[i][j] = TM[i-1][j-1]+1;
			else
				TM[i][j] = max(TM[i-1][j], TM[i][j-1]);
		}
	}
}

int main()
{
    int t, i, j;
    scanf("%d", &t);
    gets(B);

    while(t--)
    {
        gets(A);
        fill(B,B+1010,'\0');
        for(i = strlen(A)-1, j = 0; i >= 0  ;i--, j++) B[i] = A[j];
        SMC();
        cout << TM[n][m] << endl;
    }
    return 0;
}
