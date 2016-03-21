/*---------------------------------------------------------------------
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

*   O obejtivo deste código é fornecer uma implmementação
*   do AKS correta, didática e sem a utilização de
*   bibliotecas adiocinais fora do
*   padrão da língua C++

---------------------------------------------------------------------*/

/* Esta classe e suas funções são utilizadas
*  para a etapa 5 do algoritmo, a classe representa
*  um simples polinômio em uma variável arbitrária
*  e oferece as operações de atribuição e igualdade */
class polinomio{
    public:
        unsigned long long int *coef;
        unsigned int grau = 0;
        polinomio(const unsigned int &tam);
        bool operator == (const polinomio &a);
        polinomio& operator = (const polinomio &a);
        ~polinomio();
};

//Construtor da classe
polinomio::polinomio(const unsigned int &tam){
    coef = new unsigned long long int[tam+1];
    for(unsigned int i = 0; i <= tam; ++i) coef[i] = 0;
}

//Operador ==
bool polinomio::operator == (const polinomio &a){
    if(grau != a.grau) return false;
    for(unsigned int i = 0; i <= grau; ++i){ //Teste
        if(coef[i] != a.coef[i]) return false;
    }
    return true;
}

//Operador de atribuição
polinomio& polinomio::operator = (const polinomio &a){
    grau = a.grau;
    for(unsigned int i = 0; i <= grau; ++i) coef[i] = a.coef[i]; //Teste
    return *this;
}

//Destrutor
polinomio::~polinomio(){
    delete[] coef;
}

//Declaração das funções auxiliares
inline unsigned int pow_mod(const unsigned int &a, const unsigned int &b, const unsigned int &n);
inline unsigned int mdc(const unsigned int &a, const unsigned int &b);
inline unsigned int totiente_euler(unsigned int n);
inline void potencia_modular_polinomial(polinomio &resto, const unsigned int &n, const unsigned int &r, const unsigned int &a, const int log_n);

/* Teste de primalide AKS: Recebe como entrada um inteiro n positivo
*  retorna sua primalidade:  True  => Primo
*                            False => Composto
*  Para detalhes veja: < https://www.cs.auckland.ac.nz/~msta039/primality_v6.pdf > */
bool aks(unsigned int n){
    if(n < 2) return false;

    //ETAPA 1: Utilizando pessquisa binária
    unsigned int menor, maior, meio;
    double potencia;
    double logaritmo_n = log2(n);
    for(int b = 2; b <= logaritmo_n; ++b){
          menor = 1;
          maior = n;
          while(maior - menor >= 2){
                 meio = (menor + maior)/2;
                 potencia = pow((double) meio, b); //Função fornecida pelo cmath
                 if(potencia < n)  menor = meio;
                 else if(potencia > n)  maior = meio;
                 else{
                    return false;
                 }
          }
    }
    //Fim da Etapa 1

    /* Etapas 2, 3 e Etapa 4
    *  Neste passo, efetuamos a computação das 3
    *  etapas ao mesmo tempo por efeitos de otimização. */
    unsigned int r = 0, k; //'r' e a váriavel 'k' da ordem
    double logaritmo_n_2 = logaritmo_n * logaritmo_n;
    unsigned int q = 2; //Váriavel que representa iteração atual de r

    while(!r){
        if(q == n) return true; //Etapa 4

        /* Só precisamo utilizar uma chamada de mdc
        *  para as etapas 2 e 3, o valor deve ser
        *  verificado antes para etapa 4 */
        if(mdc(n, q) == 1){ //Verificando se n é relativamente primo com r atual (q)
            for(k = 1; k <= floor(logaritmo_n_2);++k){
                if(pow_mod(n, k, q) == 1) //A ordem 'k' é menor que o logaritmo
                    break;
            }
            if(k > logaritmo_n_2) //Caso contrário a ordem é algum número
                r = q;            //abaixo do logaritmo
        }
        else{ //Etapa 3
           //A iteração atual de r (q) é um fator de 'n'
            return false;
        }
        ++q;
    }
    //Fim das Etapas 2, 3 e 4

    //Etapa 5
    /*  Esta é a etapa que implementa o teorema fundamental
    *   do AKS, o método utilizado é a exponenciação binária
    *   com multiplicação quadrática. O objetivo deste código
    *   é mostrar o AKS de forma simples, qualquer operação
    *   de pré-computação para acelerar esta implementação deve
    *   ser feita aqui. Você pode tornar a implementação abaixo
    *   mais rápida pela utilização de um método rápido
    *   de multiplicação de polinômios sub-quadrática. */

    unsigned int phi; //Guarda o retorno do Totiente Euler
    phi = totiente_euler(r);

    unsigned int constante = floor(sqrt(phi)*logaritmo_n); //Limite de 'a'
    unsigned int grau_maximo = n % r; //Auxiliar para a exponenciação

    polinomio lado_dir(grau_maximo); //Construção do polinômio direito
    lado_dir.coef[grau_maximo] = 1; //P(x) = x^(n%r)
    lado_dir.grau = grau_maximo;

    polinomio lado_esq(r); //Polinomio que vai recer a computação do lado esquerdo
    for(unsigned int a = 1; a <= constante; ++a){
        lado_dir.coef[0] = a; //Polinomio direito x^(n%r) + a

        //Polinomio esq recebe o resto de ((X + a)^n / (X^r - 1)) mod n
        potencia_modular_polinomial(lado_esq, n, r, a, floor(logaritmo_n));

        //Se lado_esq <> lado_dir imprime composto
        if(!(lado_esq == lado_dir)) return false;
    }
    return true;
    //Fim Etapa 5
}

//Funções Auxiliares

//Potência modular com método binário
inline unsigned int pow_mod(const unsigned int &a, const unsigned int  &b, const unsigned int &n){
    unsigned int base, power, resp;

   resp = 1;
   base = a;
   power = b;
   while (power){
       if (power&1)
            resp = (resp*base) % n;
        power >>= 1;
        base = base%n;
        base = (base*base) % n;
   }
   return resp;
}


//Implementação do MDC de Stein
inline unsigned int mdc(const unsigned int &a, const unsigned int &b){
    unsigned int _a = a, _b = b;
    unsigned int numero_z_dir;
    if(_a == 0) return b;
    if(_b == 0) return a;
    numero_z_dir = __builtin_ctz(_a | _b); //Conta quantos 0 tem a direita do menor
    _a >>= __builtin_ctz(_a); //Desloca e atribui
    do{
        _b >>= __builtin_ctz(_b);
        if(_a > _b){
                int troca = _a;
                _a = _b;
                _b = troca;
        }
        _b = _b - _a;
    } while(_b != 0);
    return _a << numero_z_dir;
}

//Calcula a função aritmética Totiente Euler: quantidade de inteiros menores
//que a entrada e que são relativamente primos com ela
inline unsigned int totiente_euler(unsigned int n){
    unsigned long long int i;
    unsigned int resultado;
    resultado = n;

    for (i=2; i*i<=n; ++i){
        if (n % i == 0){
            while (n % i == 0){
                n /= i;
            }
            resultado -= resultado / i;
        }
    }
    if (n > 1) {
        resultado -= resultado / n;
     }
    return resultado;
}

//Pricipal cálculo da etapa 5 utilizando exponenciação binária
//e multiplicação de polinômios quadrática.
inline void potencia_modular_polinomial(polinomio &resto, const unsigned int &n, const unsigned int &r, const unsigned int &a, const int log_n){
    polinomio resultado(r); //Auxiliar, guarda resultado de muliplicacoes
    unsigned long long int novo_coeficiente;
    unsigned int anterior;
    bool bit;

    //Limpando o polinomio
    for(int i = 0; i <= r; ++i)
        resto.coef[i] = 0;

    resto.coef[0] = 1; //P(x) = 1;
    resto.grau = 0;

    for(int i = log_n; i >= 0; --i){ //Iterando sobre n binário
        //Multiplicação quadrática
        //resto = resto * resto mod (x^r - 1, n)
        resultado.grau = resto.grau * 2;
        if(resultado.grau >= r){
            resultado.grau = r - 1;
            for(unsigned int j = 0; j <= resultado.grau; ++j){
                novo_coeficiente = 0;
                for(unsigned int k = 0, l = j; k <= j; ++k, --l){ //Primeira parte
                    novo_coeficiente += resto.coef[k]*resto.coef[l];
                    novo_coeficiente %= n;
                }
                for(unsigned int k = j+1, l = r-1; k < r; ++k, --l){//Segunda parte
                    novo_coeficiente += resto.coef[k]*resto.coef[l];
                    novo_coeficiente %= n;
                }
                resultado.coef[j] = novo_coeficiente;
            }
        }
        else{ //Evita muitas multiplicações por zero
            for(unsigned int j = 0; j <= resultado.grau; ++j){
                novo_coeficiente = 0;
                for(unsigned int k = 0, l = j; k <= j; ++k, --l){ //Primeira parte
                    novo_coeficiente += resto.coef[k]*resto.coef[l];
                    novo_coeficiente %= n;
                }
                resultado.coef[j] = novo_coeficiente;
            }
        }
        resto = resultado;
        //Fim da multiplicação quadrática

        bit = (n & ( 1 << i )) >> i; //Pega o i'th bit de n
        if(bit){
            //Multiplicação por um polinômio de grau 1
            //resto = resto * (x + a) mod (x^r - 1, n)
            if(resto.grau + 1 == r) anterior = resto.coef[r - 1];
            else{
                anterior = 0;
                resto.grau += 1;
            }
            for(unsigned int j = 0; j <= resto.grau; ++j){
                novo_coeficiente = (resto.coef[j]*a) + anterior;
                anterior = resto.coef[j];
                resto.coef[j] = novo_coeficiente % n;
            }
        }
    }

    //Corrige o grau do polinômio para prox iteração
    for(unsigned int i = resto.grau; i >= 0; --i){
        if(resto.coef[i] != 0){
            resto.grau = i;
            break;
        }
    }
}
