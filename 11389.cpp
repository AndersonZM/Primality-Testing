//Anderson Zudio
//Julio Cesar
//AC 02/09/2014
#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdlib.h>
using namespace std;

int v1[101], v2[100];
int n, d, r;

int compareMA (const void * a, const void * b)
{
  return ( *(int*)a > *(int*)b );
}

int compareME (const void * a, const void * b)
{
  return ( *(int*)a < *(int*)b );
}


int main() {
	while(cin >> n >> d >> r, n){
		memset(v1,0,sizeof(v1));
		memset(v2,0,sizeof(v2));
		
		for(int i = 0; i < n; i++) 	cin >> v1[i];
		for(int i = 0; i < n; i++) 	cin >> v2[i];

  		qsort (v1, n, sizeof(int), compareMA);
  		qsort (v2, n, sizeof(int), compareME);
	
		int resp = 0;
		for(int i = 0; i < n; i++){
			int num = ((v1[i] + v2[i])-d)*r;
			resp += max(0,num);
		}
		cout << resp << endl;
	}
	return 0;
}