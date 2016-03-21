/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Testes de Primalidade Solovay-Strassen

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este código computa o teste de primalidade ultrapassado de
*   Solovay-Strassen sem suporte a inteiros de precisão arbitrária.
*   Nós recomendamos a utilização do Baillie-PSW sempre que possível.
*   O teste apresentado aqui foi utilizado para demonstrações de performance

---------------------------------------------------------------------*/

//Não aceita entradas negativas. Veja Jacobi2 no Lucas.hpp
//para um exemplo que aceita
int Jacobi(unsigned long long int a,unsigned long long int n){
    if(!a) return 0; // (0/n) = 0
    if(a==1) return 1; // (1/n) = 1
    short resp;
    int v0;
    resp=1;
    while(a){
        if(!(a&1)){
            v0 = __builtin_ctz (a);
            a>>=v0;
            if(n%8==3||n%8==5) if(v0&1) resp=-resp;
        }
        long long int aux = a;
        a = n;
        n = aux;
        if(a%4==3 && n%4==3) resp=-resp;
        a=a%n;
    }
    if(n==1) return resp;
    return 0;
}

bool solovay_strassen(const unsigned long long int & p, const int & k){
    if (!(p&1)) return false;
    unsigned long long int r, jac;
    int i;
    std::random_device rd; //Engine de randomização c++11
    std::uniform_int_distribution<unsigned long long int> uni; //distribuição uniforme em 0 ~  2^64 - 1
    std::mt19937 rgn(rd());
    for(i=0;i<k;i++){
        r = uni(rgn); //Escolhe um random
        r = 1 + (r % (p-1)); //Coloca o random no range [1, p-1]
        jac=(p+Jacobi(r,p))%p;
        if(!jac || pow_mod(r,(p-1)/2,p)!=jac){
            return false;
        }
    }
    return true;
}
