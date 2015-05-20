//Anderson zudio
//Bruno Masquio
//Victor Cracel
//AC em 26/08/2014

#include<iostream>

using namespace std;

typedef struct comb{
	int a, b, c, s;
}comb;

comb poss[100];

bool verifica(int a1, int a2, int a3, int b1, int b2, int b3, int c1, int c2, int c3){
	return a1==b1 || a1==b2 || a1==b3 || a2==b1 || a2==b2 || a2==b3 || a3==b1 || a3==b2 || a3==b3 ||
	       a1==c1 || a1==c2 || a1==c3 || a2==c1 || a2==c2 || a2==c3 || a3==c1 || a3==c2 || a3==c3 ||
	       c1==b1 || c1==b2 || c1==b3 || c2==b1 || c2==b2 || c2==b3 || c3==b1 || c3==b2 || c3==b3;
}

int main(){
	int n, caso = 0;
	while(cin >> n, n){
		int resposta = -1;
		for(int i = 0; i < n; i++)
			cin >> poss[i].a >> poss[i].b >> poss[i].c >> poss[i].s;
			
		for(int i = 0; i < n; i++){
			for(int j = i+1; j < n; j++){
				for(int k = j+1; k < n; k++){
					int a1 = poss[i].a, a2 = poss[i].b, a3 = poss[i].c,
						b1 = poss[j].a, b2 = poss[j].b, b3 = poss[j].c,
						c1 = poss[k].a, c2 = poss[k].b, c3 = poss[k].c;
						if(!verifica(a1, a2, a3, b1, b2, b3, c1, c2, c3)) 
							resposta = max(resposta, poss[i].s + poss[j].s + poss[k].s);
				}
			}
		}

		cout << "Case " << ++caso << ": " << resposta << endl;
	}
	
	return 0;
}