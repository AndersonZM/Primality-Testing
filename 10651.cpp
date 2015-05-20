//Anderson Zudio
//Victor Cracel
//Bruno Masquio
//AC em 26/08/2014

#include <iostream>
#include <string>
#define MAX 50;

using namespace std;

int resposta;

void jogo(string meu_jogo){
	int contador = 0;
	
	auto jogo_original = meu_jogo;
	
	for(int i = 0; i < 11; i++){
		if(meu_jogo[i] == 'o' && meu_jogo[i + 1] == 'o'){
			if(i + 2 < 12 && meu_jogo[i + 2] == '-'){
				meu_jogo[i + 2] = 'o';
				meu_jogo[i] = meu_jogo[i + 1] = '-';
				jogo(meu_jogo);
				meu_jogo = jogo_original;
			}
			
			if(i - 1 >= 0 && meu_jogo[i - 1] == '-'){
				meu_jogo[i - 1] = 'o';
				meu_jogo[i] = meu_jogo[i + 1] = '-';
				jogo(meu_jogo);
				meu_jogo = jogo_original;
			}
		}
	}
	
	for(int i = 0; i < 12; i++)
		if(meu_jogo[i] == 'o') contador++;
		
	resposta = min(contador,resposta); 
}

int main() {
	int n;
	cin >> n;
	while(n--){
		string linha;
		cin >> linha;
		
		resposta = MAX;
		
		jogo(linha);
		
		cout << resposta << endl;
		
	}
	return 0;
}