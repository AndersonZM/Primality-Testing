/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Testes de Primalidade Miller-Rabin

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este código computa o teste de primalidade Miller-Rabin muito
*   utilizado na prática. Não há suporte para inteiros de
*   precisão arbitrária. O teste de Baillie-PSW e melhor para todas
*   as entradas possíveis que este código computa, a menos que você execute.
*    este código com somente 2 iterações (baixa probabilidade de primo).

---------------------------------------------------------------------*/

bool miller_rabin(const unsigned long long int & p, const int &k) {
    if ((p&1)){ //
        unsigned long long int p1 = p-1;
        unsigned long long int v0 = __builtin_ctz (p1), p2; // número de 0 a direita
        unsigned long long int r,s,mod;
        std::random_device rd; //Engine de randomização c++11
        std::uniform_int_distribution<unsigned long long int> uni; //distribuição uniforme em 0 ~  2^64 - 1
        std::mt19937 rgn(rd());

        // removando potencias de 2 de p1
        p2 = p1 >> v0; // removi todos os 0 a direita
        for(int i=0;i < k; ++i){
            r = uni(rgn); //Escolhe um random
            r = 1 + (r % (p-1)); //Coloca o random no range [1, p-1]
            s = p2;
            mod = pow_mod (r,s,p); // mod = (r^s) mod p
            while (s!=(p1) && mod!=1 && mod!=(p1)){
                // voltou para a condição original ou
                //encontrou uma das duas condições de paradas do item
                mod = pow_mod (mod,2,p);
                s *= 2;
            }
            if (mod!=p1 && !(s&1)){
                return false;
            }
        }
    }
    else {return false;}
    return true;
}
