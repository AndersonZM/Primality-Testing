/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Potência modular com exponenciação binária

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este código é utilizado para a execução dos testes de primaidade
*   probabilísticos.
---------------------------------------------------------------------*/

inline unsigned long long int pow_mod(const unsigned long long int &a, const unsigned long long int  &b, const unsigned long long int &n){
    __int128 base, power, resp, aux; //utilize long long para sistemas que nao suportam 128 bits

   resp = 1;
   base = a;
   power = b;
   while (power){
       if (power&1){

            aux = (resp%n)*(base%n);
            resp = aux % n;
       }
        power >>= 1;
        base = base%n;
        aux = base*base;
        base = aux % n;
   }

   return resp;
}
