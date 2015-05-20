//AC em 26/04/2014
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;
long int TM[1001][1001], A[110], B[110], N1, N2;

void SMC(){
    int i, j;
    for (i=1; i<=N1; i++) TM[i][0]=0;
    for (j=1; j<=N2; j++) TM[0][j]=0;
	for (i=1; i<=N1; i++)
		for (j=1; j<=N2; j++)
			if (A[i-1] == B[j-1]) TM[i][j] = TM[i-1][j-1]+1;
			else  	              TM[i][j] = max(TM[i-1][j], TM[i][j-1]);
}

int main()
{
    int cont = 0;
    while(cin >> N1, cin >> N2)
    {
        if (N1 == 0 && N2 == 0) break;
        cont++;
        for(int k = 0; k < N1; k++)cin >> A[k];
        for(int j = 0; j < N2; j++)cin >> B[j];
    SMC();
    cout <<"Twin Towers #"<< cont << endl << "Number of Tiles : "<<TM[N1][N2] << endl << endl;

    }

    return 0;
}
