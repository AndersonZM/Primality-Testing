/*--------------------------------------------------------------------
    Universidade do Estado do Rio de Janeiro
    Autores: Anderson Zudio de Moraes
             Victor Cracel Messner
*   Testes de Primalidade de Lucas

*   O código apresentado faz parte do projeto de conclusão de curso
*   dos autores pelo título de Bacharel, apresentado ao
*   Intituto de Matemática e Estatística,
*   da Universidade do Estado do Rio de Janeiro.
*   Você pode modificar e distribuir este código livremente
*   através dos termos do GNU General Public como publicado
*   pelo Free Software Foundation através da versão 3 da licença,
*   ou qualquer versão subsequente por opção sua.

*   Este código computa o teste de primalidade de Lucas para o
*   Baillie-PSW. O teste é inútil sozinho e não deve ser utilizado
*   na prática. Há muitas entradas que vão apontar como primos
*   sendo estas compostas, elas são chamadas de pseudoprimos de lucas.

---------------------------------------------------------------------*/


//Esta função aceita entradas negativas
inline int Jacobi2(long long int a,long long int n){
    if(!a) return 0;
    if(a==1) return 1;
    short resp;
    int v0;
    resp=1;
     if(a<0){
        a=-a;
        if(n%4==3) resp=-resp;
    }
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

bool lucas (long long int p){
    __int128 d,ud,P,Q,U,V,U2,V2,aux,uaux,n;
    short sinal,jac;
    ud=5; sinal = 1;
    while (1){
        d = ud * sinal;
        jac = Jacobi2(d,p);
        if (!jac && ud!=p){
            return false;
        }
        if (jac == -1){
            break;}
        sinal = -sinal;
        ud+=2;
    }
    P=1;
    Q= (1-d)/4;
    U=0; V=2;
    U2=1; V2=P;
    n = (p+1)/2;
    while (n){
        U2 = (p+(U2*V2))%p;
        V2 = (p+(V2 * V2 - Q - Q))%p;
        if (n&1){
            uaux = U;
            U = U2 * V + U * V2;
            if ((U&1))
            U+=p;
            U/=2;
            V = V2*V + U2 * uaux * d;
            if ((V&1))
            V+=p;
            V/=2;
            if (U<0) U-=((U/p)+1)*p;
            if (V<0) V-=((V/p)+1)*p;
            U%=p;
            V%=p;
        }
        Q = (p + (Q * Q)) % p;
        n/=2;
    }
    if (U) return false;
    return true;
}
