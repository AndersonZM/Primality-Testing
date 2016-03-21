/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Sistema de prova de primalidade - ECPP

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este é um exemplo de uso do sistema ECPP. O código segue o algoritmo
*   ECPP de Atkin. A biblioteca utilizada  o PARI, vista em
*   < http://pari.math.u-bordeaux.fr/doc.html >
*   Aumente a restrição de memória n MAIN para testes grandes.
*   Este é possivelmente o melhor sistema para porvar a primalidade de um número
*   conhecido atualmente. Não esqueça de utilizar -std=c++11 e opcionalmente
*   -Ofast na compilação.

---------------------------------------------------------------------*/
#include <iostream> //I/O
#include <pari/pari.h> //biblioteca PARI
#include <cmath> //Operações de matematica
#include <random> //Engine de random do c++11
#include <chrono>   //Contagem de tempo
#include <getopt.h> //opções
#include <random> //Random engine do c++

/* Contém informações básicas de parametros de calibramento
*  Você deve calibrar o DMAX, a pre-computação e o número de
*  pontos a serem testados. Também a precisão é um argumento a ser calibrado
*  recomendamos no mínimo 1000 casas decimais. */
#include "small_primes.h"

bool certificado;

//Retorna um fator provavelmente primo de m maior que limite ou NULO
GEN provavelmente_fatorado(GEN m, GEN limite);

/* Algoritmo que computa o polinomio de classe Hilbert,
   Veja See "A Course in Computational Algebraic Number Theory"
   de Henri Cohen pagina 415. Dado um discriminante negativo D
   o algoritmo computa o polinômio de grau h(D) em Zn[X]
   tal que j((D+sqrt(D))/2) é uma raiz. */
 GEN somatorio(GEN q);
 GEN hilbert(long D);


/* Verifica se tem como encontrar a curva elíptica  E(a, b)
   pela definição de Lenstra, sabendo que #E(a, b) = m. */
bool curva_ideal(long D, GEN n, GEN *raiz, long *qtd_raiz, GEN *a, GEN *b);

/* Computa P' = P*k na curva elítpica E(a, b) pela definição de Lenstra */
int computa_P(std::pair<GEN, GEN> P, GEN k, GEN n_i, GEN a, std::pair<GEN, GEN> *Pr);

/* Não há validação de entrada em n.
   A função aguarda um intiero positivo de
   precisão arbitrária.
   Se o sistema retornar FALSE, há uma pequena chance do sistema
   ter desistido de provar a primalidade do número. */
bool ecpp(const std::string n){
    pari_sp av = avma, av2, av3, av4;   //Ponteiros para garbage collection
    GEN n_i = gp_read_str(n.c_str()), q, m, a, b, raiz, invariante, c, g, x, y, n_i1, teste;
    std::pair<GEN, GEN> P, P1, P2;
    GEN U, V; //Recebe as solucoes de cornacchia
    long D, qtd_raiz, k, i = 0;
    bool achei;

    while(gcmpgs(n_i, SMALL_PRIME) == 1){
        av2 = avma;

        /* Se a execução for um número abaixo de 2^64 - 1
           então o Baillie-psw é determinístico e você pode parar de computar.
           Retorne o resultado de Baillie-PSW.
           Este código não para aqui pois o objetivo é demonstrar o sistema ECPP*/
        if(!BPSW_psp(n_i)){ //Provavelmente primo
                av = avma; //Garbage collection
                return false;
        }

        for(D = -3; D >= -D_MAX; D--){
            achei = false;

            //Testando discriminante fundamental
            if(D%4 != 0 && D%4 != -3 && D%4 != 1) continue;

            av3 = avma; //Garbage collection


            GEN D_MOD_Ni = addsi(D, n_i);
            //Computa o simbolo de Jacobi, a entrada  provavelmente primo
            if(kronecker(D_MOD_Ni, n_i) != 1){
                    avma = av3;
                    continue;
            }

            //Retorna a solucao da eq diogfantina para as variaveis A e B
            GEN GEN_d = stoi(D);
            if(!cornacchia2(negi(GEN_d), n_i, &U, &V)){
                    avma = av3;
                    continue;
            }

            //Determinando m com fator q bom
            teste = itor(n_i, PRECISION);
            teste = sqrtnr(teste, 4); teste = addrs(teste, 1);
            teste = powru(teste, 2);

            //m = n_i + 1 + U
            n_i1 = addii(n_i, gen_1); //Auxiliar que guarda N_i + 1
            m = gcopy(n_i1);
            m = addii(m, U);
            q = provavelmente_fatorado(m, teste);

            if(q == NULL){
                //m = n_i + 1 - U
                affii(n_i1, m); m = subii(m, U);
                q = provavelmente_fatorado(m, teste);
            }

            //Especiais
            if(q == NULL){
                if(D == -4){
                    //m = n_i + 1 + 2*V
                    affii(n_i1, m); m = addii(m, mulis(V, 2));
                    q = provavelmente_fatorado(m, teste);

                    if(q == NULL){
                        //m = n_i + 1 - 2*V
                        affii(n_i1, m); m = subii(m, mulis(V, 2));
                        q = provavelmente_fatorado(m, teste);
                   }

                }
                else if(D == -3){
                    //m = n_i + 1 + (U + 3*V)/2
                    affii(n_i1, m); m = addii(m, divii(addii(U, mulis(V, 3)),gen_2));
                    q = provavelmente_fatorado(m, teste);

                   if(q == NULL){
                        //m = n_i + 1 - (U + 3*V)/2
                        affii(n_i1, m); m = subii(m, divii(addii(U, mulis(V, 3)),gen_2));
                        q = provavelmente_fatorado(m, teste);
                    }

                    if(q == NULL){
                        //m = n_i + 1 + (U - 3*V)/2
                        affii(n_i1, m); m = addii(m, divii(subii(U, mulis(V, 3)),gen_2));
                        q = provavelmente_fatorado(m, teste);
                    }

                    if(q == NULL){
                        //m = n_i + 1 + (U - 3*V)/2
                        affii(n_i1, m); m = subii(m, divii(subii(U, mulis(V, 3)),gen_2));
                        q = provavelmente_fatorado(m, teste);
                    }
                }
            }

            //Nao consegui fatorar m
            if(q == NULL){
                    avma = av3; //Garbage Collection
                    continue;
            }

            //Eu tenho o #E(a, b), vou tentar construir a curva
            qtd_raiz = 0;
            for(int tipo = 0; tipo < 2; tipo++){
                if(!tipo){ //Casos especiais
                    if (D == -3) {
                        a = Fp_red(gen_0, n_i);
                        b = Fp_red(gen_m1, n_i);
                        qtd_raiz = 1;
                    }
                    else if (D == -4) {
                        a = Fp_red(gen_m1, n_i);
                        b = Fp_red(gen_0, n_i);
                        qtd_raiz = 1;
                    }
                }
                else if(!curva_ideal(D, n_i, &raiz, &qtd_raiz, &a, &b)) continue;

                for(int tentativa = 0; tentativa < qtd_raiz; tentativa++){
                    if(tipo){
                        invariante = Fp_red(gel(raiz, tentativa+1), n_i);
                        c = Fp_red( mulii(invariante, ginvmod(subis(invariante, 1728), n_i)), n_i);
                        a = Fp_red(mulis(c, -3), n_i);
                        b = Fp_red(mulis(c, 2), n_i);
                    }

                    //Consegui achar uma curva, vou computar g para a multiplicação complexa
                    while(true){
                        do{
                            g = randomi(n_i);
                        }
                        while(gequal0(g));
                        if(kronecker(g, n_i) == -1){

                            if(D == -3) //Caso especial
                                if(gequal1(Fp_pow(g, diviuexact(subiu(n_i, 1), 3), n_i)))
                                    continue;
                                break;
                        }
                    }

                    //Vou tentar encontrar um ponto para satisfazer o teorema
                    for(int t_ponto = 0; t_ponto < PTRIES; t_ponto++){

                        //Encontrando um ponto aleatório
                        y = gen_0;
                        while(gequal0(y)){
                            do {
                                do
                                    x = randomi(n_i);
                                while (gequal0(x));
                                y = powiu(x, 3);
                                y = addii(y, mulii(x, a));
                                y = addii(y, b);
                                y = Fp_red(y, n_i);
							} while (kronecker(y, n_i) == -1);
							y = Fp_sqrt(y, n_i);
							if(y == NULL) y = gen_0;
						}

                        P.first = gcopy(x); P.second = gcopy(y);
						k = 0;
                        while(true){
                            //computando P1 e P2
                            int respP1 = computa_P(P, m, n_i, a, &P1);
                            if(respP1 == 1){
                                //Nao foi possivel realizar a computaçao
                                //n_i é composto, desistindo
                                //de provar a primalidade de n_0
                                avma = av;
                                return false;
                            }

                            int respP2 = computa_P(P, divii(m, q), n_i, a, &P2);
                            if(respP2 == 1){
                                //Nao foi possivel realizar a computaçao
                                //n_i é composto, desistindo
                                //de provar a primalidade de n_0
                                avma = av;
                                return false;
                            }

                            if(respP1 == -1 && respP2 == 0){
                                //Achei um ponto que atende o teorema
                                achei = true;
                                break;
                            }
                            ++k;

                            //Vou tentar denovo com outra curva
                            if(D == -3){
                                if(k >= 6) break;
                                b = gmul(b, g);
                            }
                            else if(D == -4){
                                if(k >= 4) break;
                                a = gmul(a, g);
                            }
                            else{
                                if(k >= 2) break;
                                a = gmul(gmul(g, g), a);
                                b = gmul(gpowgs(g, 3), b);
                            }
                            a = Fp_red(a, n_i);
                            b = Fp_red(b, n_i);
                        }
                        if(achei) break;
                    }
                    if(achei) break;
                }
                if(achei) break;
            }
            if(achei) break;
            else avma = av3; //Garbage collection
        }

        if(-D == D_MAX+1){
                //Tentei muito, melhor desistir
                //Não quer dizer que n_0 é composto,
                //eu estou aqui porque cheguei no meu limite
                //e não consigo provar a primalidade de n_i
                avma = av;
                return false;
        }

        //Imprimidno o certificado
        if(certificado){
            std::cout << "n_" << i;
            pari_printf(" = %Ps\n", n_i);
            std::cout << "E(a, b), a = ";
            pari_printf("%Ps, b = %Ps\n", a, b);
            std::cout << "m = ";
            pari_printf("%Ps, q = %Ps\n", m, q);
            std::cout << "P = (";
            pari_printf("%Ps, %Ps)\n", P.first, P.second);
            std::cout << "P1 = (";
            pari_printf("%Ps, %Ps)\n", P1.first, P1.second);
            std::cout << "P2 = (";
            pari_printf("%Ps, %Ps)\n", P2.first, P2.second);
        }

        i++;
        n_i = gerepileupto(av2, gcopy(q));  //Garbage Collection
    }

    //n_I é pequeno, vou fazer uma pesquisa binária nos primos que conheço
    long n_long = itos(n_i);

    long low = 0, high = LAST_PRIME, mid;
    while(low <= high){
        mid = (low + high)/2;
        if(primos[mid] == n_long){
                av = avma;
                if(certificado)
                    std::cout << n_long << " e um primo determinado pelo AKS." << std::endl;
                return true; //Consegui provar que é primo
        }
        else if(primos[mid] < n_long) low = mid + 1;
        else high = mid - 1;
    }

    //Restaurando
    av = avma;
    return false; //Pela natureza do baillie-psw, eu não devo executar isto
}

GEN provavelmente_fatorado(GEN m, GEN limite){
    pari_sp ltop = avma;
    GEN q = gcopy(m);
    GEN gen_3 = addii(gen_2, gen_1), gen_5 = addii(gen_3, gen_2);

    while(gequal0(modii(q, gen_2))) q = diviiexact(q, gen_2);
    while(gequal0(modii(q, gen_3))) q = diviiexact(q, gen_3);
    while(gequal0(modii(q, gen_5))) q = diviiexact(q, gen_5);


    if(gcmp(q, limite) != 1 || gequal(q, m)){
            avma = ltop;
            return NULL;
    }


    if(BPSW_psp(q)) return gerepileupto(ltop, q);

    avma = ltop;
    return NULL;
}

GEN somatorio(GEN q){
    pari_sp av = avma;
    GEN erro = gen_1;
	long e1, e2, n = 1, sinal = -1;
	GEN termo1, termo2, soma = gen_0;

    for (long i = 0; i < 1000; i++)
		erro = gdivgs(erro, 2);

	do {
		e1 = n * (3 * n - 1) / 2;
		e2 = n * (3 * n + 1) / 2;
		termo1 = gpowgs(q, e1);
		termo2 = gpowgs(q, e2);
		termo1 = gadd(termo1, termo2);
		termo1 = gmulgs(termo1, sinal);
		soma = gadd(soma, termo1);
		sinal = -sinal;
		n++;
	} while (gcmp(gabs(termo1, PRECISION), erro) == 1);

	return gerepileupto(av, gadd(gen_1, soma));
}


GEN hilbert(long D){
    pari_sp av = avma;
    GEN gen_D, c_poly, q_poly, ca, cb, tau, J, c, pi2, q, q2, deltaq, deltaq2, f, f1, JJ;
    long b, B, t, a;
    gen_D = gsqrt(stor(D, PRECISION), PRECISION);

    c_poly = cgetg(3, t_POL);
    c_poly[1] = evalvarn(0);
    gel(c_poly, 2) = gen_1; //as contas acabam gerando t_COMPLEX
    c_poly = normalizepol(c_poly);

    b = D % 2;
    if(b < 0) b += 2;
	B = (long) sqrt(labs(D) / 3.0);

    do{
        t = (b * b - D) / 4;
        a = (b > 1) ? b : 1;

        do{
            if (t % a == 0){
                ca = mulss(2, a); ca = itor(ca, PRECISION);
                cb = stor(-b, PRECISION);
                tau = gdiv(gadd(cb, gen_D), ca);
                pi2 = constpi(PRECISION);
                pi2 = mulrs(pi2, 2);
                c = cgetc(PRECISION);
                gel(c, 1) = itor(gen_0, PRECISION);
                gel(c, 2) = pi2;

                q = gexp(gmul(c, tau), PRECISION);
                q2 = gmul(q, q);

                deltaq = gmul(q, gpowgs(somatorio(q), 24));
                deltaq2 = gmul(q2, gpowgs(somatorio(q2), 24));

                f = gdiv(deltaq2, deltaq);
                f1 = gmulsg(256, f);
                f1 = gadd(gen_1, f1);
                f1 = gpowgs(f1, 3);
                J = gdiv(f1, f);

                if (a == b || a * a == t || b == 0)
                    q_poly = mkpoln(2, gen_1, gneg(J));
                else{
                    JJ = gadd(gmul(greal(J), greal(J)), gmul(gimag(J), gimag(J)));
                    q_poly = mkpoln(3, gen_1, gmulsg(-2, greal(J)), JJ);
                }

                c_poly = RgX_mul(c_poly, q_poly);
            }
            a++;
        } while(a * a <= t);

        b += 2;
    } while(b <= B);

    //O codigo abaixo não é muityo bom, mas resolve vários bugs
    for(long i = 2; i <= degree(c_poly)+2; i++){
        if(typ(gel(c_poly, i)) == 6) gel(c_poly, i) = ground(greal(gel(c_poly, i)));
        else gel(c_poly, i) = ground(gel(c_poly, i));
    }

    return gerepileupto(av, gcopy(c_poly));
}

bool curva_ideal(long D, GEN n, GEN *raiz, long *qtd_raiz, GEN *a, GEN *b){
    GEN hpoly, menor_coeficiente;
    *qtd_raiz = 0;
    pari_sp av = avma;
    hpoly = hilbert(D);

    menor_coeficiente = gel(hpoly, 2);
    if(!Z_ispower(menor_coeficiente, 3)){
        avma = av;
        return false;
    }
    hpoly = RgX_to_FpX(hpoly, n);
    *raiz = gerepileupto(av, FpX_roots(hpoly, n));
    *qtd_raiz = (lg(*raiz)-1);

    //Se eu nao tenho raízes, deve retornar falso
    if(*qtd_raiz > 0) return true;
    else return false;
}

int computa_P(std::pair<GEN, GEN> P, GEN k, GEN n_i, GEN a, std::pair<GEN, GEN> *Pr){
    pari_sp ltop = avma, lbot;
    int resposta = 1; //-1 se infinito, 0 se nao tem divisor e 1 se tem
    std::pair<GEN, GEN> A, O, C;
    GEN d, d_x, d_y, m, inv;


    A.first = gcopy(P.first); A.second = gcopy(P.second);
    O.first = gen_0; O.second = gen_1;

    while(resposta && (cmpis(k, 0) == 1)){
        if(!gequal0(Fp_red(k, gen_2))){
            //Fast Point mult & addition
            d = gcdii(Fp_red(subii(O.first, A.first), n_i), n_i);
            k = subis(k, 1);
            if(gequal1(d) || gequal(d, n_i)) resposta = 1;
            else resposta = 0;
            if(!(gequal0(A.first) && gequal1(A.second))){
                if(gequal0(O.first) && gequal1(O.second)){
                    O.first = gcopy(A.first); O.second = gcopy(A.second);
                }
                else if(resposta){
                    d_x = Fp_red(subii(O.first, A.first) , n_i);
                    d_y = Fp_red(subii(O.second, A.second), n_i);
                    if(gequal(A.first, O.first) && gequal0(Fp_red(addii(A.second, O.second) , n_i))){
                        C.first = gen_0; C.second = gen_1;
                    }
                    else{
                        inv = Fp_invsafe(d_x, n_i);
                        if(inv == NULL) inv = gen_0;
                        m = Fp_red(mulii(d_y, inv), n_i);
                        C.first = Fp_red(subii(mulii(m, m), addii(A.first, O.first)), n_i);
                        C.second = Fp_red(subii(mulii(m, subii(A.first, C.first)), A.second), n_i);
                    }
                    O.first = gcopy(C.first); O.second = gcopy(C.second);
                }
            }


        }
        else{
            d = gcdii(Fp_red(mulis(A.second, 2), n_i), n_i);
            k = diviuexact(k, 2);
            if(gequal1(d) || gequal(d, n_i)) resposta = 1;
            else resposta = 0;

            if(resposta){
                inv = Fp_invsafe(mulsi(2, A.second), n_i);
                if(inv == NULL) inv = gen_0;
                m = Fp_red(mulii(addii(mulis(mulii(A.first, A.first), 3), a), inv), n_i);
                C.first = Fp_red(subii(mulii(m, m), mulis(A.first, 2)), n_i);
                C.second = Fp_red(subii(mulii(m, subii(A.first, C.first)), A.second), n_i);
                A.first = gcopy(C.first); A.second = gcopy(C.second);
            }
        }
    }

    lbot = avma;
    Pr->first = Fp_red(O.first, n_i);
    Pr->second = gerepile(ltop, lbot, Fp_red(O.second, n_i));

    if(gequal0(Pr->first) && gequal1(Pr->second)) return -1;
    return !resposta;
}

//Declarações de variáveis para a função main
std::chrono::steady_clock::time_point inicio, fim;
std::chrono::duration<double> tempo_decorrido;
std::string n; //Entrada
bool teste, tempo, individual, encontra;
int prox_opcao, k;
const char* const opcoes = "ietc";
const char* nome_programa;

void inline configura(int, char**);
void inline uso();

const struct option argumentos[] = {
	{"individual", 0, NULL, 'i'},
	{"encontra", 0, NULL, 'e'},
	{"tempo", 0, NULL, 't'},
	{"certificado", 0, NULL, 'c'},
	{NULL, 0,NULL, 0}
};

int main(int argc, char **argv){
    nome_programa = argv[0];
    if(argc == 1) uso();
    else configura(argc, argv);

    pari_init(100000000, 65536); //Aumente para entradas maiores

    /*--------------------------------------
    Tenta provar a primalidade da entrada
    --------------------------------------*/
    if(individual && !encontra){
        while(std::cin >> n){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = ecpp(n);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << "Provei que e primo " << n << ". ";
            else std::cout <<"Nao consegui provar " << n << ". ";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
    }

    /*-------------------------------------------------
    Me dê o número de bits, eu te dou um inteiro primo
    --------------------------------------------------*/
    else if(!individual && encontra){
        while(std::cin >> k){
            std::random_device rd;
            std::uniform_int_distribution<short> uni;
            std::mt19937 rgn(rd());
            pari_sp av = avma;
            GEN primo = gen_0, gen_R;
            gen_R = powis(gen_2, k-1);

            inicio = std::chrono::steady_clock::now();
            do{
                av = avma;
                setrand(stoi(uni(rgn)));
                primo = addii(gen_R, randomi(gen_R));
                if(!mpodd(primo)) continue;
                n = std::string(itostr(primo));
                if(ecpp(n)) break;
                avma = av;
            }while(true);
            fim = std::chrono::steady_clock::now();
            pari_printf("%Ps", primo);
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
    }
    else uso();

    pari_close();
    return 0;
}


/*------------------------------------------
O código abaixo é somente para configuração
-------------------------------------------*/

void inline uso(){
    std::cout << "Utilize assim: " << nome_programa << " opcao [-t] [-c] [ < arquivo_entrada ]" << std::endl;
    std::cout <<
			"	-i	--individual	Prova a primalidade de primos individuais\n"
			"	-e	--encontra  	Encontra um primo de n bits.\n"
			"	-t	--tempo			Ativa a contagem de tempo.\n"
			"	-c	--certificado   Imprime o certificado.\n" ;
	exit(0);
}

void inline configura(int argc, char** argv){
    certificado = individual = tempo = encontra = false;
    do {
		prox_opcao = getopt_long(argc, argv, opcoes, argumentos, NULL);
		switch (prox_opcao)
		{
			case 'i':
				individual = true;
				break;
			case 'c':
				certificado = true;
				break;
			case 'e':
				encontra = true;
				break;
			case 't':
				tempo = true;
				break;
			case -1:
				break;

			default:
				uso();
		}
	}
	while (prox_opcao != -1);
}


