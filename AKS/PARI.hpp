/*------------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Teste de Primalidade AKS - PARI

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
*   utiliza a biblioteca PARI, que pode ser vista
*   em < http://pari.math.u-bordeaux.fr/doc.html >

-----------------------------------------------------------------------*/

#include <pari/pari.h> //código da biblioteca principal

GEN potencia_modular_polinomial(const GEN &base, const GEN &mod, const long &poly_mod, const ulong logaritmo_n); //Computa potencia polinomial via algoritmo binário

/* Observe que a função abaixo recebe a entrada
*  em forma de STRING. Ela está esperando um string
*  que contenha somente um número inteiro natural,
*  válido para a execução do AKS. Não há validação
*  de entrada no código. */
bool aks_pari(const std::string _n){ //Evita o programdor a utilizar I/O do PARI
    /*----------------------------------------------------------*
    | Nao esqueca de iniciar a STACK PARI para as funcoes       |
    | PARI trabalha com uma stack propria que exige que         |
    | o programador faca garbage collection por conta propria   |
    | o ponteiro abaixo serve para retornar a stack para        |
    | onde ela estava depois da execucao de alguma funcao       |
    *----------------------------------------------------------*/
    pari_sp av = avma, av2; //Guarda o ponteiro da stack, recupera depois
    GEN n; //variável do tipo long long int *
           //implmenta TODOS os tipos internos do PARI
    n = gp_read_str(_n.c_str()); //Não faz validação de entrada

    //Etapa 1
    if(Z_isanypower(n, NULL)){
        avma = av; //Garbage collection, restaura o topo da stack
        return false;
    }
    //Fim Etapa 1

    //Etapas 2, 3 e 4
	bool teste;
    long k; //WARNING! Possivel overflow para n > 2^(32768*raiz(2))
	GEN r = gen_2;
	GEN n_REAL = itor(n, 3); //Variavel que guarda o valor de n como t_REAL "para o log"
	double logaritmo_n =  dbllog2r(n_REAL); //Precia ser real, se nao for causa erros de calculo
    double logaritmo_n_2 = logaritmo_n*logaritmo_n;

    while(1){
        if(gequal(r, n)){  //Etapa 4
            avma = av; //Garbage collection, restaura o topo da stack
            return true;
        }
        if(gequal1(gcdii(n, r))){ //Verificando se n é relativamente primo com q
            teste = true;
            for(k = 1; k <= logaritmo_n_2; k++){
                if(gequal1(Fp_powu(n, k, r))){
                    teste = false;
                    break;
                }
            }
            if(teste) break; //Ordem menor que logaritmo_n_2
        }
        else{ //Etapa 3
            avma = av; //Restaura o ponteiro
            return false;

        }
        r = gaddgs(r, 1); //r++
    }
    //Fim da Etapa 2, Etapa 3 e Etapa 4*/

    //Etapa 5
	//-- Cuidado: N muito grande pode causar overflow no grau do polinomio --
	//-- Limitacao da biblioteca! Grau max (2^31-1)  ... r vale no maximo log2(N)^5 ...
    //WARNING! Possivel Overflow para n > 2^(raiz[5](2^32 - 1)) (~23 dígitos)
	long r_long = gtolong(r);

    //Limite dessas variaveis dependem de (r), se r nao estiver em Overflow as contas sao seguras
    long phi; //Funcao totiente de Euler
    phi = eulerphiu(r_long);
    long constante = floor(sqrt(phi)*logaritmo_n); //piso de raiz(phi(n)) log2 (n)
    GEN n_mod_r = Fp_red(n, r);
    long n_mod_r_long = gtolong(n_mod_r);
    GEN lado_direito, lado_esquerdo;

    //Aqui temos o polinomio do lado direito sendo construido
    //com funções de baixo nível da biblioteca, cuidado
    //ao alterar o código abaixo. Possível point exception
    lado_direito = cgetg(n_mod_r_long + 3, t_POL);
    lado_direito[1] = evalvarn(0);
    for(ulong x = 2; x <= n_mod_r_long+2; ++x) gel(lado_direito, x) = gen_0;
    gel(lado_direito, n_mod_r_long + 2) = gen_1; //x^(r%n)
    lado_direito = normalizepol(lado_direito);

    for(long a; a <= constante; ++a){
        av2 = avma; //Guarda o ponteiro da stack
        lado_esquerdo = mkpoln(2, gen_1, gen_0); //cria o lado esquerdo, safe
        gel(lado_esquerdo, 2) = stoi(a); //x + a

        //Ceil falha só quando n é uma potencia de 2
        lado_esquerdo = potencia_modular_polinomial(lado_esquerdo, n, r_long, ceil(logaritmo_n));
        gel(lado_esquerdo, 2) = subis(gel(lado_esquerdo, 2), a); //lado_esquerdo -= a;

        //O código abaixo é uma verificação em baixo nível equivalente a
        //lado_direto != lado_esquerdo ? composto : continua;
        if(degree(lado_esquerdo) != degree(lado_direito)){
            avma = av; //Garbage collection
            return false;
        }
        for(ulong x = 0; x <= degree(lado_esquerdo); ++x){
            if(!gequal(gel(lado_esquerdo, x+2), gel(lado_direito, x+2))){
                avma = av; //Garbage collection
                return false;
            }
        }

        avma = av2; //Garbage collection, restaura o topo da stack
    }
    avma = av; //Garbage collection, restaura o topo da stack
    return true;
    //Fim Etapa 5
}

/* Esta função é uma tradução da função de mesmo nome do AKS simples.
*  Nós somente traduzimos o código em C para código em PARI.
*  O PARI é extremamente eficiente aqui, ele utiliza uma
*  uma multiplicação de polinômios sub-quadrática. */
GEN potencia_modular_polinomial(const GEN &base, const GEN &mod, const long &poly_mod, const ulong logaritmo_n){
    GEN resto = pol_1(0);
    GEN pow2i, ndivpow, truncate, bit;

    for(ulong i = logaritmo_n; i > 0; i--){
        resto = FpX_sqr(resto, mod); //resto = resto*resto mod n
        resto =  ZX_mod_Xnm1(resto, poly_mod); //resto em Zn[x]/x^r - 1

        //Pegando o i-ésimo bit de n
        pow2i = powuu(2, i-1);
        ndivpow = rdivii(mod, pow2i, 3);
        truncate = floorr(ndivpow);
        bit = Fp_red(truncate, gen_2);

        if(gequal1(bit)){
            resto = FpX_mul(resto, base, mod); //resto = resto*(x + a)
            resto = ZX_mod_Xnm1(resto, poly_mod); //resto em Zn[x]/x^r - 1
        }

    }
    resto = RgX_to_FpX(resto, mod); //Reduz o polinomio para Zn[x]
    return resto;
}
