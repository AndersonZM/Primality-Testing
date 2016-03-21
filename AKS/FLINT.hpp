/*------------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Teste de Primalidade AKS - FLINT

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
*   utiliza a biblioteca FLINT, que pode ser vista
*   em < http://www.flintlib.org/ >

-----------------------------------------------------------------------*/
#include <flintxx.h>
#include <fmpzxx.h>
#include <fmpz_mod_polyxx.h>
/* -------------------------------------------
--Header Default do FLINT inclui por default:
--<stdlib.h>, <stdio.h>, <flint.h>, <gmp.h>
------------------------------------------- */


using namespace flint;


/* Dessa vez a função recebe um fmpz
*  que é o inteiro de precisão arbitrária do flint,
*   o input é simples com interface de c++ do flintxx:
*   cin >> n;                                           */
bool aks_flint(const fmpz &n){
    //Etapa 1
    mpz_t mpz_n; // tipo utilizado pelo GMP, utilizado na primeira etapa
    mpz_init(mpz_n);
    fmpz_get_mpz(mpz_n, &n);
    //O FLINT nao oferece uma função para
    //computar a etapa 1. O GMP oferece:
    if(mpz_perfect_power_p(mpz_n)){
        return false;
    }

	//Fim da Etapa 1

	//Etapas 2, 3 e 4
	bool teste;
    long k; //Overflow para n > 2^(32768*raiz(2))
	fmpz r = 2, pot_mod, gcd; //N muito grande pode causar overflow no r na etapa 5
	double logaritmo_n =  fmpz_dlog(&n)/fmpz_dlog(&r); //Precia ser real, se nao for causa erros de calculo
    double logaritmo_n_2 = logaritmo_n*logaritmo_n;

    while(1){
        if(r == n){  //Etapa 4
            return true;
        }
        fmpz_gcd(&gcd, &n, &r);
        if(fmpz_is_one(&gcd)){ //Verificando se n é relativamente primo com q
            teste = true;
            for(k = 1; k <= logaritmo_n_2; ++k){
                fmpz_powm_ui(&pot_mod, &n, k, &r);
                if(fmpz_is_one(&pot_mod)){
                    teste = false;
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
	//-- Cuidado: N muito grande pode causar overflow no grau do polinomio  fmpz_mod_poly_t --
	//-- Limitacao da biblioteca! Grau max (2^31-1) ... r vale no maximo log2(N)^5 ...
	//WARNING! Possivel Overflow para n > 2^(raiz[5](2^32 - 1)) (~23 dígitos)
    fmpz phi, coeficiente; //Funcao totiente de Euler
    fmpz_euler_phi(&phi, &r);
    double phi_d = fmpz_get_d(&phi);
    long constante = floor(sqrt(phi_d)*logaritmo_n); //piso de raiz(phi(n)) log2 (n)

    //Aqui nós montamos o polinomio do lado dir e esq com as ferramentas do flint
    slong r_long = fmpz_get_si(&r);
    fmpz_mod_poly_t lado_direito, lado_esquerdo, modulo;
    fmpz_mod_poly_init(lado_direito, &n);
    fmpz_mod_poly_init(lado_esquerdo, &n);

    //Aqui nós montamos o polinomio do modulo
    fmpz_mod_poly_init(modulo, &n);
    fmpz n_1 = n - 1;
    fmpz_mod_poly_set_coeff_ui(modulo, r_long, 1);
    fmpz_mod_poly_set_coeff_fmpz(modulo, 0, &n_1); //p(x) = x^r + (n-1)
                                                   //computar ocm este p(x) é equivalente
                                                   //a computar com p'(x) = x^r - 1
    fmpz n_mod_r = n % r;
    slong n_mod_r_long = fmpz_get_si(&n_mod_r);
    fmpz_mod_poly_set_coeff_ui(lado_direito, fmpz_get_si(&n_mod_r), 1); //x^(n%r)

    for(slong a; a <= constante; ++a){
        fmpz_mod_poly_set_coeff_ui(lado_esquerdo, 1, 1);
        fmpz_mod_poly_set_coeff_ui(lado_esquerdo, 0, a);
        //A função abaixo computa o p(x) = (x + a)^n mod (x^r + n-1, n)
        fmpz_mod_poly_powmod_fmpz_binexp(lado_esquerdo, lado_esquerdo, &n, modulo);

        fmpz_mod_poly_get_coeff_fmpz(&coeficiente, lado_esquerdo, 0);
        coeficiente -= a;
		fmpz_mod_poly_set_coeff_fmpz(lado_esquerdo, 0,  &coeficiente); //lado_esquerdo -= a

		//if(lado_esquerdo != lado_direito)
		if(!fmpz_mod_poly_equal(lado_esquerdo, lado_direito)) return false;
         fmpz_mod_poly_realloc(lado_esquerdo, 0); //Clear lado_esquerdo
    }
    return true;
	//Fim Etapa 5
}
