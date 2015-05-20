//Anderson Zudio
//AC em 31/07/2014
#include<iostream>
#include<cmath>
#include <vector>
#include <unordered_map>


using namespace std;

unordered_map <string, int> pessoas;
vector <string> ordem;
bool flag = false;

int main()
{
    int n, gift, qtd, pagou;
    string nome;
    while(cin >> n){
        pessoas.clear();
        ordem.clear();
        for(int i = 0; i < n; i++){
            cin >> nome;
            pessoas.insert({nome, 0});
            ordem.push_back(nome);
        }
        for(int i = 0; i < n; i++){
            cin >> nome >> gift >> qtd;
            pagou = 0;
            for(int j = 0; j < qtd; j++){
                string p_nome;
                cin >> p_nome;
                pessoas[p_nome] += gift/qtd;
                pagou += gift/qtd;
            }
            pessoas[nome] -= pagou;
        }

        if(flag) cout << endl;
        else flag = true;
        for(auto i = ordem.begin(); i != ordem.end(); i++)
            cout << *i << " " << pessoas[*i] << endl;
    }

	return 0;
}
