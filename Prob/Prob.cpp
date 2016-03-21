/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Testes de Primalidade Monte Carlo - Miller-Rabin, Solovay-Strassen
*   e Baillie-PSW
*   OBS: Para as entradas possíveis desse código, o Baillie-PSW é
*   determinístico.

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este é um exemplo de função principal para executar
*   os testes de primalidade probabilsticos apresentados no projeto.
*   Não esqueça que este código utiliza a biblioteca Chrono para
*   a contagem de tempo. Compile com a opção -std=c++11,
*   recomendamos o uso de -Ofast (opcional).

---------------------------------------------------------------------*/

#include <iostream> //cin; cout
#include <ratio> //duration
#include <chrono> //contagem de tempo
#include <getopt.h> //opções
#include <string> //Algumas funções para lidar com string
#include <cmath> //Utilizado pelas funções de logaritmo, exponenciação e raiz
#include <random> //Random engine do c++

#include "POW_MOD.hpp" //pow_mod(), utilizado por todos
#include "SS.hpp" //Solovay-Strassen
#include "MR.hpp" //Miller-Rabin
#include "Lucas.hpp" //Lucas do Baillie-PSW
#include "BPSW.hpp" //Baillie-PSW

std::chrono::steady_clock::time_point inicio, fim;
std::chrono::duration<double> tempo_decorrido;
unsigned long long int n; //Entrada
bool teste, SS, MR, BPSW, tempo, quiet;
int prox_opcao, k;
const char* const opcoes = "smbtq";
const char* nome_programa;

void inline configura(int, char**);
void inline uso();

const struct option argumentos[] = {
	{"solovay-strassen", 0, NULL, 's'},
	{"miller-rabin", 0, NULL, 'm'},
	{"baillie-psw", 0, NULL, 'b'},
	{"tempo", 0, NULL, 't'},
	{"quiet", 0, NULL, 'q'},
	{NULL, 0,NULL, 0}
};

int main(int argc, char* argv[]){
    nome_programa = argv[0];
    if(argc == 1) uso();
    else configura(argc, argv);

    /*-------------------
      Solovay-Strassen
    -------------------*/
    if(SS && !MR && !BPSW){
        if(!quiet){
            std::cout << "Solovay-Strassen -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

        while(std::cin >> n >> k){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = solovay_strassen(n, k);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n << " e provavelmente primo.";
            else std::cout << n << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
    }

    /*----------------
        Miller-Rabin
    ----------------*/
    else if(!SS && MR && !BPSW){
        if(!quiet){
            std::cout << "Miller-Rabin -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

	    while(std::cin >> n >> k){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = miller_rabin(n, k);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n << " e provavelmente primo.";
            else std::cout << n << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
	    }
    }

    /*--------------
       Baillie-PSW
    ---------------*/
   else if(!SS && !MR && BPSW){
        if(!quiet){
            std::cout << "Baillie-PSW -- Contagem de Tempo: ";
            tempo ? std::cout << "ativada.\n" : std::cout << "desativada.\n";
        }

        while(std::cin >> n){
            tempo_decorrido = std::chrono::steady_clock::duration::zero();
            inicio = std::chrono::steady_clock::now();
            teste = baillie_psw(n);
            fim = std::chrono::steady_clock::now();
            if(teste) std::cout << n << " e primo."; //deterministico até 2^64
            else std::cout << n << " e composto.";
            if(tempo){
                tempo_decorrido = std::chrono::duration_cast<std::chrono::duration<double>>(fim - inicio);
                std::cout << " Tempo de execucao: " << tempo_decorrido.count();
            }
            std::cout << std::endl;
        }
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
			"	-s	--solovay-strassen	Utiliza o teste Solovay-Strassen\n"
			"	-m	--miller-rabin 		Utiliza o teste Miller-Rabin.\n"
			"	-b	--baillie-psw	    Utiliza o teste Baillie-PSW.\n"
			"	-t	--tempo			Ativa a contagem de tempo.\n"
			"	-q	--quiet			So imprime o resultado, minimalista.\n" ;
	exit(0);
}

void inline configura(int argc, char** argv){
    quiet = tempo = MR = SS = BPSW = false;
    do {
		prox_opcao = getopt_long(argc, argv, opcoes, argumentos, NULL);
		switch (prox_opcao)
		{
			case 'm':
				MR = true;
				break;
			case 's':
				SS = true;
				break;
			case 'b':
				BPSW = true;
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


