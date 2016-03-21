/*------------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Teste de Primalidade AKS - NTL

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   O objetivo deste código é fornecer uma implmementação
*   do AKS correta, rápida e que suporte entradas
*   maiores que a implementação simples. Este código
*   utiliza a biblioteca NTL, que pode ser vista
*   em < http://www.shoup.net/ntl/ >

-----------------------------------------------------------------------*/

#include <gmp.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pX.h>
/* -------------------------------------------
--Header Default do NTL inclui por default:
--<cstdlib>, <cmath> e <iostream>
------------------------------------------- */

inline void totiente_euler(long &resultado, long n);

/* Observe que a função abaixo recebe a entrada
*  em forma de STRING. Ela está esperando um string
*  que contenha somente um número inteiro natural,
*  válido para a execução do AKS. Não há validação
*  de entrada no código. */
bool aks_ntl(const std::string n){ //Evita o programdor a utilizar I/O do NTL
	//Etapa 1
    mpz_t mpz_n; // tipo utilizado pelo GMP
    mpz_init(mpz_n);
    mpz_set_str(mpz_n, n.c_str(), 10);
    //O NTL nao oferece uma função para
    //computar a etapa 1. O GMP oferece:
    if(mpz_perfect_power_p(mpz_n))
            return false;
	//Fim da Etapa 1

	//Etapas 2, 3 e 4
	NTL::ZZ ZZ_n;  // tipo utilizado pelo NTL, inteiro prec arbitrária
	ZZ_n = to_ZZ(NTL::conv<NTL::ZZ>(n.c_str())); //Não há validação de entrada
	bool teste;
    long k; //Overflow para n > 2^(32768*raiz(2))
	NTL::ZZ r(2), n_mod_r; //N muito grande pode causar overflow no r na etapa 5

	double logaritmo_n = log(ZZ_n)/log(2); //Precisa ser real, se nao for causa erros de calculo
    double logaritmo_n_2 = logaritmo_n*logaritmo_n;

    while(1){
        if(r == ZZ_n){  //Etapa 4
            return true;
        }
        if(IsOne(GCD(ZZ_n, r))){ //Verificando se n é relativamente primo com q
            teste = true;
            rem(n_mod_r, ZZ_n, r);
            for(k = 1; k <= logaritmo_n_2; ++k){
                if(IsOne(PowerMod(n_mod_r, k, r))){ // (n^k)%r==1 --> n%r passado para evitar BAD_ARGUMENTS
                    teste = false;                  // resolve um problema da biblioteca
                    break;
                }
            }
            if(teste) break; //Ordem menor que logaritmo_n_2
        }
        else //Etapa 3
            return false;
        ++r;
    }
    //Fim da Etapa 2, Etapa 3 e Etapa 4

	//Etapa 5
	//-- Cuidado: N muito grande pode causar overflow no grau do polinomio ZZ_pX --
	//-- Limitacao da biblioteca! Grau max (2^31-1) ... r vale no maximo log2(N)^5 ...
	//WARNING! Possivel Overflow para n > 2^(raiz[5](2^32 - 1)) (~23 dígitos)
    long phi; //Funcao totiente de Eulerĵ
    long r_long = NTL::conv<long>(r);
    totiente_euler(phi, r_long);
    long int constante = floor(sqrt(phi)*logaritmo_n); //piso de raiz(phi(n)) log2 (n)

    NTL::ZZ_p::init(ZZ_n);//Ativa a modularidade com n para o NTL
    NTL::ZZ_pX f(r_long, 1); f -= 1; //x^r - 1, polinomio para acelerar as contas
    const NTL::ZZ_pXModulus modulo(f);  //internas do NTL com pre-computacao
    NTL::ZZ_pX lado_direito(NTL::conv<long>(n_mod_r), 1); //p(x) = x^(n mod r)
    NTL::ZZ_pX lado_esquerdo(1, 1); //p(x) = x

    for(long a; a <= constante; ++a){
        SetCoeff(lado_esquerdo, 1); lado_esquerdo += a;//p(x) = x + a
        PowerMod(lado_esquerdo, lado_esquerdo, ZZ_n, modulo); //(lhs = (x + a)^n mod (x^r - 1, n)
        lado_esquerdo -= a;

        if(lado_esquerdo != lado_direito) return false;
        clear(lado_esquerdo);
    }
    return true;
	//Fim Etapa 5

}

//Esta função é exatamente a mesma utilizada
//no AKS simples, ela é bem rápida
inline void totiente_euler(long &resultado, long n){
    resultado = n;

    for (long i = 2; i*i<=n; ++i){
        if ((n % i) == 0){
            while ((n % i) == 0){
                n /= i;
            }
            resultado -= resultado / i;
        }
    }
    if (n > 1) {
        resultado -= resultado / n;
     }
}
