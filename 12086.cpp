//Anderson Zudio
//Raquel Marcolino
//AC em 9/5/2015

#include <iostream>
#include <vector>
#include <cstdio>
#include<string.h>
//Uma fórmula de debug útil que pode ser utilizada
#define RSQPrint(A,B) printf("RSQ(%d,%d) = %d\n",A,B,RSQ(A,B))
using namespace std;

int N, pot[20010], x, y;
char c[5];
vector<int> V,FT;

//Função LSOne(), mais explicações sobre ela podem ser
//encontradas no slide.
int LSOne(int x){
    return (x)&(-x);
}

//A função abaixo calcula o RSQ para um resultado individual
//olhando os índices de sua representação em binário
//Exemplo: 6 na base 10 = 110 na base 2, ou seja,
//RSQ(6) = FT[6] + FT[4]
int RSQ(int b){
    int soma=0;
    while(b){
        soma+=FT[b];
        b-=LSOne(b);
    }
    return soma;
}

//A função abaixo calcula o RSQ para um dado intervalo
//Se um RSQ individual de a dá como resposta o RSQ(1,a),
//então um RSQ(1,b) - RSQ(1,a-1) deve dar a resposta de
//RSQ(a,b). Observe que a está incluído no intervalo.
int RSQ(int a, int b){
    if(a==1){
        return RSQ(b);
    }else{
        return RSQ(b) - RSQ(a-1);
    }
}

//Essa função ajusta o valor de um determinado índice k para
//qualquer valor v. Ela ajusta também quaisquer índices que
//estejam sob a responsabilidade deste.
//Exemplo: 5 na base 10 = 101 na base 2, ou seja,
//ajustar o índice 5 ajusta também o índice 6 (110) e o
//índice 8 (1000), mas isso depende também do tamanho da árvore.
void Ajusta(int k, int v){
    int tam;
    tam=FT.size();
    while(k<tam){
        FT[k]+=v;
        k+=LSOne(k);
    }
}

int main() {
    int t,i,cont, casei = 0, valores[200100];
    bool first = true;
    cont= 0; //Essa variável conta a quantidade de números do input.
    /** Input recomendado:
        2 4 5 5 6 6 6 7 7 8 9
    **/
    while(scanf("%d", &N) == 1){
        FT.assign(N+1,0);
        if(N == 0) break;
        if(!first) printf("\n");
        else first = false;
        memset(valores, 0, sizeof(valores));
        for(auto i = 1; i <= N; i++){
            auto a = 0;
            scanf("%d", &a);
            valores[i] = a;
            Ajusta(i, a);
        }


        /*for(i=0;i<cont;++i){
            //Ajustando os valores na árvore de acordo com o vetor inicial.
            Ajusta(V[i],1);
        }*/

        printf("Case %d:\n", ++casei);

        while(scanf("%s", c)){
            if(c[0] == 'E') break;
            scanf("%d %d", &x, &y);
            if(c[0] ==  'M'){
                printf("%d\n", RSQ(x, y));
            }
            else{
                auto temp = y - valores[x];
                Ajusta(x, temp);
                valores[x] = y;
            }
        }
    }

    /*while(cin >> t){
        //Inserindo elementos no vetor inicial.
        V.push_back(t);
        ++cont;
    }
    //Vamos iniciar um vetor com índices de 1 a 10 para a Árvore de Fenwick.
    //Lembrando que 1 a 10 é a faixa dos valores do vetor inicial.
    FT.assign(11,0);
    for(i=0;i<cont;++i){
        //Ajustando os valores na árvore de acordo com o vetor inicial.
        Ajusta(V[i],1);
    }
    cout << "RSQ(1,1) = " << RSQ(1,1) << endl; //FT[1] = 0
    Ajusta(5,2); //Somou 2 aos índices 5, 6 e 8
    cout << "RSQ(1,10) = " << RSQ(1,10) << endl; //FT[10] + FT[8] = 1 + 12 = 13
    return 0;*/
}
