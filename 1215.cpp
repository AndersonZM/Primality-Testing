//Anderson Zudio
//AC em 29/03/2015

#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace std;

int t, c, resp,  letras[10001][26];
string s;
vector<int> cortes(1001);
set<int> range;

int main() {
    cin >> t;
    for(int i = 0; i < t; i++) {
        cortes.clear();
        cin >> c;
        for(int i = 0; i < c; i++){
            int corte_atual;
            cin >> corte_atual;
            cortes.push_back(corte_atual);
        }

        cin >> s;
        for(int i = 0; i <= s.length(); i++){
            for(int j = 'a'; j <= 'z'; j++){
                if(!i) letras[i][j] = 0;
                else{
                    if(s[i-1] == j)
                        letras[i][j-'a'] = letras[i-1][j-'a'] + 1;
                    else
                        letras[i][j-'a'] = letras[i-1][j-'a'];
                }
            }
        }

       /* for(int i = 0; i <= s.length(); i++){
            for(int j = 0; j < 26; j++)
                cout << letras[i][j] << " ";
            cout << endl;
        }*/

        range.clear();
        range.insert(0);
        range.insert(s.length());
        resp = 0;
        for(int i = 0;  i < c; i++) {
            auto corte_atual = cortes.at(i);
            auto it = range.lower_bound(corte_atual);
            auto corte1 = *it;
            auto corte2 = *--it;
            range.insert(corte_atual);

           /* for(auto it2 = cortes.begin(); it2 != cortes.end(); it2++)
                cout << *it2 << " ";
            cout << endl;*/

            for(int j = 'a'; j <= 'z'; j++) {
                auto lado_esquerdo = letras[corte_atual][j-'a'] - letras[corte2][j-'a'];
                auto lado_direito = letras[corte1][j-'a'] - letras[corte_atual][j-'a'];
                if (lado_esquerdo > 0 ^ lado_direito > 0) resp++;
            }
        }

        cout << resp << endl;
    }
    return 0;
}
