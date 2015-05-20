//Anderson Zudio Moraes
//Julio Cesar
//Victor Cracel Messner
//AC 02/05/2014

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>

/* Este problema foi resolvido com mudanças no algoritmo de SCSM.
A idéia é a mesma para quando a subseq chega em valores menor que zero.
Se o valor da subseq é maior que zero, verifique se ele é maior que X.
Se ele for maior que X, comece a retirar todos os valores da subseq a partir do primeiro termo afim de achar alguma subseq
DENTRO da subseq que você acabou de computar, como há negativos na entrada,
pode ser que exista uma subseq dentro da subseq que tenha o valor maior que
X e seja menor do que a subseq que você acabou de achar.

*/

using namespace std;

int S[500005], n, x, tamanho;

void SCSM(){
    long long int ss, sst;
    int atamanho, atamanhot, cs, cst;
	ss = atamanho = cs = 0;

	for( int i = 1; i <= n; i++ ) {

        //o inteiro ss guarda a soma atual, se ss + próximo termo for maior que zero
        //atualiza o tamanho atual, esta idéia é exatamente a mesma coisa do SCSM original
		if( ( ss + S[i] ) >= 0 ) {
			ss += S[i];
			atamanho++;
			if( cs == 0 ) cs=i;
			if( ss >= x ) {
			    //Caso a minha subseq for maior que X, verifique se retirando os primeiros termos eu consigo achar uma
			    //subseq dentro da subseq atual que é menor. Isto é possível devido aos negativos ou em caso de S[i] ser muito grande tal que seja perto de X.
                    sst = ss;
                    atamanhot = atamanho;
                    cst = cs;
                    while(cst <= i)
                    {
                        //Atualizando o tamanho
                        if(sst >= x)
                        {
                            if(tamanho == -1) tamanho = atamanhot; //Se o código nunca entrar aqui, a seq toda era negativa. Imprimir -1
                            else tamanho = min(tamanho, atamanhot);
                        }
                        sst -= S[cst++];
                        atamanhot--;
                    }
            }
        }
        //Mesma idéia original SCSM
		else {
			ss = cs = atamanho = 0;
        }
	}
}

int main()
{
    int t;
    cin >> t;
    while(t--)
    {
        cin >> n >> x;
        for(int i = 1; i <= n; i++)
            cin >> S[i];

        tamanho = -1;

            SCSM();

        cout << tamanho << endl;
    }
}
