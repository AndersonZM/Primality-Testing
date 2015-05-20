//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
//AC em 28/04/2014

/* Problema Área:
    Time Limit: 1 sec
    Solução O(n^4) não é viavel, é necessária uma solução com complexidade menor que isso ou Time Limit.
    O problema foi resolvido computando uma matriz com a soma em linha (pode ser implementado em coluna também) de todos os preços.
    Com esta nova matriz temos como computar o preço de uma submatriz de maneira mais eficiente, para calcular o preço de qualquer submatriz basta somar os valores na última coluna:
    Por exemplo, o preço para a11 até aij onde i = 3 e j = 3 pode ser obtido de a13 + a23 + a33.
    Caso eu queira o preço de uma submatriz dentro da matriz acima, basta diminuir os valores das primeiras linhas:
    Por exemplo, o preço para a12 até aij onde i =3 e j= 3 pode ser obtida de a23 + a33. Utilizando a área computada no exemplo acima: (a13 + a23 + a33) - a13
                                                                                                                                        área anterior - a13
                                                                                                                                        basta computar a subtração de a13.
    Com isto a solução fica fácil: basta subtrair o valor da soma até achar algo menor ou igual ao orçamento da entrada e guardar a maior área.


    Exemplo prático:
       Seja a entrada  1 2 3      matriz soma em linha:       1  3  6
                       4 5 6                                  4  9  15
                       7 8 9                                  7  15 24
                       Preço de a11 até a33: 6 + 15 + 24 = 45              para a22 até a33:  15 + 24  (anterior menos a primeira linha)
                                             (1+2+3 + 4+5+6 + 7+8+9)

*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>

using namespace std;

int main()
{
    //Declarações
    int t, cont, S, P, n, m, K, i, j, k, l, Satual, Patual, g, h, M[110][110];
    cin >> t;
    cont = 0;
    while(t--)
    {

        cin >> n >> m >> K;

        for(i = 1; i <= n; i++)
            for(j = 1; j <= m; j++)
                cin >> M[i][j];


        //Guardando a matriz soma por linha.
        for(i = 1; i <= n; i++)
            for(j = 2; j <= m; j++)
                 M[i][j] += M[i][j-1];


        //Computando resposta:
        S = P = 0;
        for(j = 1; j <= m; j++)
        {

            //Para cada loop, retirar os valores da primeira coluna
            if(j > 1)
                for(g = 1; g <= n; g++)
                    for(h = j; h <= m; h++)
                        M[g][h] -= M[g][j-1];
            //Após este loop, a coluna j atual foi eliminada do problema e já foi computada

            //O inteiro k identifica a coluna do começo da submatriz que está sendo computada
            for(k = j; k <= m; k++)
            {
                Satual = Patual = 0;

                //O inteiro i identifica a linha do começo da submatriz a ser computada
                for(i = 1, l = 1; i <= n; i++)
                {
                    //Adiciona o preço que está armazena em aik, este preço contém o valor da linha toda até a coluna k
                    Patual += M[i][k];
                    //Se este preço for maior do que meu orçamento, subtraia das primeiras linhas até encaixar no meu orçamento
                    if(Patual > K)
                    {
                        while(Patual > K)
                        {
                            Patual -= M[l++][k];
                        }
                    }

                    //Calculando a área da minha submatriz que começa em alj e termina em aik, atualiza a resposta se necessário
                    Satual = (i - l + 1)*(k - j + 1);
                    if(Satual > S)
                    {
                        S = Satual;
                        P = Patual;
                    }
                    else if(Satual == S)
                        if(Patual < P) P = Patual;
                }
            }

        }


        cout << "Case #" << ++cont << ": " << S << " " << P << endl;
    }
    return 0;
}
