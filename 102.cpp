//Anderson Zudio
//Victor Cracel
//AC em 27/08/2014
#include <iostream>
using namespace std;

int main() {
	int b1, g1, c1, b2, g2, c2, b3, g3, c3, movimento, movimento_atual;
	string resposta;
	
	while(cin >> b1 >> g1 >> c1 >> b2 >> g2 >> c2 >> b3 >> g3 >> c3){
		//GCB
		movimento_atual = b1 + c1 + g2 + b2 + g3 + c3;
		movimento = movimento_atual; 
		resposta = "GCB";
		
		//GBC
		movimento_atual = b1 + c1 + g2 + c2 + b3 + g3;
		if(movimento_atual <= movimento){
			movimento = movimento_atual;
			resposta = "GBC";
		}
		
		//CGB
		movimento_atual = b1 + g1 + b2 + c2 + g3 + c3;
		if(movimento_atual <= movimento){
			movimento = movimento_atual;
			resposta = "CGB";
		}
		
		//CBG
		movimento_atual = b1 + g1 + g2 + c2 + b3 + c3;
		if(movimento_atual <= movimento){
			movimento = movimento_atual;
			resposta = "CBG";
		}
		
		//BGC
		movimento_atual = g1 + c1 + b2 + c2 + b3 + g3;
		if(movimento_atual <= movimento){
			movimento = movimento_atual;
			resposta = "BGC";
		}
		
		//BCG
		movimento_atual = g1 + c1 + b2 + g2 + b3 + c3;
		if(movimento_atual <= movimento){
			movimento = movimento_atual;
			resposta = "BCG";
		}
		
		
		cout << resposta << ' ' << movimento << endl;
	}
	return 0;
}