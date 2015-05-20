//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 11/04/2014

#include <iostream>
#include <string>
#include <vector>

using namespace std;

struct no
{
    int x, mod, pai;
};

int n, resto[100010], digitos[10], achei, a, cont;
string f;
no aux;

int main()
{
    vector<no> vector1;
    vector<no>::iterator vectori;
    while(cin >> f >> n)
    {
		if(!n) break;
        vector1.clear();
        fill(resto,resto+100009,0);
        fill(digitos,digitos+10,0);
        achei = 0;

        for(int i = 0; i < f.length(); i++)
            digitos[f[i] - '0'] = 1;

		for(int i = 1; i < 10; i++)
        {
            if(digitos[i])
            {
                aux.x = i;
                aux.pai = 0;
                if(i % n) aux.mod = i % n;
                else
                {
                    achei = i;
                    break;
                }
                resto[aux.mod] = 1;
                vector1.push_back(aux);
            }
        }
        if(achei)
        {
            cout << achei << endl;
            continue;
        }

        vectori = vector1.begin();
        cont = 1;
        while((achei != 1) && (vectori != vector1.end()))
        {
            a = vectori->mod*10;
            for(int i = 0; i < 10; i++)
            {
                if(digitos[i])
                {
                    aux.x = i;
                    aux.pai = cont;
                    aux.mod = (a + aux.x) % n;
                    if(!resto[aux.mod])
                    {
                        vector1.push_back(aux);
                        resto[aux.mod] = 1;
                        if(!aux.mod)
                        {
                            achei = 1;
                            break;
                        }
                    }
                }
            }
            vectori = vector1.begin() + cont;
            cont++;
        }

        /*vectori = vector1.begin();
        cont = 1;
        while(vectori != vector1.end())
        {
            cout << vectori->x << " " << vectori->mod << " " << vectori->pai << " " << cont << endl;
            cont++;
            vectori++;
        }
        cout << endl;*/

        if(achei)
        {
            f = "";
            vectori = vector1.end();
            vectori--;
            while(vectori->pai != 0)
            {
                f += vectori->x + '0';
                vectori = vector1.begin() + vectori->pai - 1;
            }
            f += vectori->x + '0';
            for(int i = f.size()-1; i >= 0; i--)
                cout << f[i];
            cout << endl;

        }
        else cout << "impossible" << endl;

    }
    return 0;
}
