/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Teste de Primalidade AKS - Simples

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este é um exemplo de função principal para executar
*   as versões do AKS apresentadas no projeto. Não esqueça
*   que este código utiliza a biblioteca Chrono para a contagem
*   de tempo. Compile com a opção -std=c++11, recomendamos o uso
*   de -Ofast (opcional).

---------------------------------------------------------------------*/

#include <iostream> //cin; cout
#include <ratio> //duration
#include <chrono> //contagem de tempo
#include <getopt.h> //opções
#include <string> //Algumas funções para lidar com string
#include <cmath> //Utilizado pelas funções de logaritmo, exponenciação e raiz

#include "AKS.hpp" //AKS simples
#include "NTL.hpp" //AKS NTL
#include "FLINT.hpp" //AKS FLINT
#include "PARI.hpp"  //AKS PARI

std::chrono::steady_clock::time_point inicio, fim;
std::chrono::duration<double> tempo_decorrido;
unsigned int n_int; //Entrada do AKS simples
std::string n_string; //Entrada do NTL e do PARI
fmpz n_fmpz; //Entrada do FLINT
bool teste, simples, ntl, pari, flin, tempo, quiet;
int prox_opcao;
const char* const opcoes = "spnftq";
const char* nome_programa;

void inline configura(int, char**);
void inline uso();

const struct option argumentos[] = {
	{"simples", 0, NULL, 's'},
	{"pari", 0, NULL, 'p'},
	{"ntl", 0, NULL, 'n'},
	{"flint", 0, NULL, 'f'},
	{"tempo", 0, NULL, 't'},
	{"quiet", 0, NULL, 'q'},
	{NULL, 0,NULL, 0}
};

int main(int argc, char* argv[]){
    nome_programa = argv[0];
    if(argc == 1) uso();
    else configura(argc, argv);

    /*-------------
      AKS SIMPLES
    --------------*/
    if(simples && !pari && !ntl && !flin){
        if(!quiet){
            std::cout << "AKS Simples -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

        while(std::cin >> n_int){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = aks(n_int);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n_int << " e primo.";
            else std::cout << n_int << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
    }

    /*-------------
        AKS PARI
    --------------*/
    else if(!simples && pari && !ntl && !flin){
        if(!quiet){
            std::cout << "AKS PARI -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }
        /*-------------------------------------------------
	     Inicia o stack para pre computacoes do PARI.
         Primeiro argumento eh o numero de bytes que o pari
         tem para trabalhar -- 5mb utilizado,
         o segundo eh o primo maximo que ele pode
         pre-computar -- 2^16 utilizado
	    -------------------------------------------------*/
	    pari_init(100000000, 65536);

	    while(std::cin >> n_string){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = aks_pari(n_string);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n_string << " e primo.";
            else std::cout << n_string << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
	    }
        //Retorna a memoria utilizada pelo PARI STACK para o sistema operacional
        pari_close();
    }

    /*-------------
        AKS NTL
    --------------*/
   else if(!simples && !pari && ntl && !flin){
        if(!quiet){
            std::cout << "AKS NTL -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

        while(std::cin >> n_string){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = aks_ntl(n_string);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n_string << " e primo.";
            else std::cout << n_string << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
    }

    /*-------------
        AKS FLINT
    --------------*/
    else if(!simples && !pari && !ntl && flin){
        if(!quiet){
            std::cout << "AKS FLINT -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

        while(std::cin >> n_fmpz){ //interface do flint para c++
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = aks_flint(n_fmpz);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n_fmpz << " e primo.";
            else std::cout << n_fmpz << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
	    flint_cleanup(); //Libera a memória utilizada pelo FLINT
    }
    else uso();

    return 0;
}

/*------------------------------------------
O código abaixo é somente para configuração
-------------------------------------------*/

void inline uso(){
    std::cout << "Utilize assim: " << nome_programa << " opcao [-t] [-q] [ < arquivo_entrada ]" << std::endl;
    std::cout <<
			"	-s	--simples		Utiliza o AKS Simples\n"
			"	-p	--pari  		Utiliza o AKS com PARI.\n"
			"	-n	--ntl			Utiliza o AKS com NTL.\n"
			"	-f	--flint 		Utiliza o AKS com FLINT.\n"
			"	-t	--tempo			Ativa a contagem de tempo.\n"
			"	-q	--quiet			So imprime o resultado, minimalista.\n" ;
	exit(0);
}

void inline configura(int argc, char** argv){
    quiet = tempo = ntl = simples = flin = pari = false;
    do {
		prox_opcao = getopt_long(argc, argv, opcoes, argumentos, NULL);
		switch (prox_opcao)
		{
			case 's':
				simples = true;
				break;
			case 'p':
				pari = true;
				break;

			case 'n':
				ntl = true;
				break;
            case 'f':
				flin = true;
				break;
			case 't':
				tempo = true;
				break;
			case 'q':
				quiet = true;
				break;
			case -1:
				break;

			default:
				uso();
		}
	}
	while (prox_opcao != -1);
}
