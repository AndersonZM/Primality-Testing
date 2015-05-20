//Anderson Zudio
//Bruno Masquio
//Giancarlo França
//Raquel Marcolino
//AC em 05/08/2014
//C++11

#include <iostream>
#include <set>
#include <cmath>

#define BIT 64


using namespace std;

set<unsigned long long, less<unsigned long long> > myset; //Set com ordenação começando pelo menor número

/*Para resolver este problema basta achar todos os números que tenham como base um inteiro elevado a um outro inteiro não primo.
Nesta implementação foi usado um set, o set é um container que permite que a inserção seja ordenada e os elementos sejam únicos.*/


//Pow da biblioteca cmath do C++11 retorna double e tem como entrada a base sendo double e o expoente também.
//Como estamos lidando com números muito grandes (extamente o tamanho de unsigned l l), a melhor opção foi usar um pow próprio
unsigned long long ipow( unsigned long long base, unsigned long long expo )
{
    unsigned long long oi, resp=1;
    for(oi = 0; oi < expo; oi++)
            resp*=base;

    return resp;
}

bool TestaPrimo (int q)
{    int rq;   bool primo=true;
     rq = sqrt(q);
     for (int i=2; i<=rq; i++)
         if ((q % i) == 0)
         {   primo = false;
             break;
         }
     return (primo);
}


int main() {

	unsigned long long i,j,power, powerbase;
	powerbase=2;
	myset.insert(1); //Set permite inserção ordenada e com elementos únicos

	while(powerbase < 65537){ //Aqui o loop foi superestimado, pouca perda de desempenho pois para i = 2 o valor de power chega no máximo a 63. i = 3, power = 40 e para cada incremento em i o valor de power desce muito até chegar em 0.
		power = ceil(BIT  / (log(powerbase)/log(2)));//Aqui é calculado o expoente máximo
		for(i = power; i>0; i--){
			if(!TestaPrimo(i)){ //Teste de primalidade, se o expoente é primo, não tentar inserir no set

				unsigned long long x = ipow(powerbase, i);

				if(64*log(2) > i*log(powerbase)) myset.insert(x); //Para casos que estouram 2^64-1
			}
		}
		powerbase++;
	}

	//Set já está ordenado pela função less<unsigned long long>
	for(auto k = myset.begin(); k!= myset.end(); k++){
		cout << *k << endl;
	}

	return 0;
}
