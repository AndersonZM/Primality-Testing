//Anderson Zudio
//AC em 30/07/2014
#include <iostream>
using namespace std;

int main() {
	int juice, mile, t, n, i, j, aux, nn;
	cin >> t;
	for(i = 1; i <= t; i++){
		cin >> n;
		mile = juice = 0;
		for(j = 0; j < n; j++){
			cin >> nn;
			aux = nn;
			while(aux >= 0){
				aux -= 30;
				mile++;
			}
			aux = nn;
			while(aux >= 0){
				aux -= 60;
				juice++;
			}
		}
		mile *= 10;
		juice *= 15;
		if(juice > mile) cout << "Case " << i << ": Mile " << mile << endl;
		else if(mile > juice) cout << "Case " << i << ": Juice " << juice << endl;
		else cout << "Case " << i << ": Mile Juice " << mile << endl;
	}
	return 0;
}
