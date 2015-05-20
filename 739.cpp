//AC em 22/01/2013
//Anderson Zudio de Moraes
#include <string>
#include <iostream>

using namespace std;

string T, Soundex;

string geraSoundex(string T)
{
    int n, m , cont;
    n = T.size();
    cont = 0;
    string S;
    S = T[0];
    m = S.size();
    for(int i = 1; i < n; i++)
    {
        if (m == 4) break;
        if (((T[i] == 'B') || (T[i] == 'P') || (T[i] == 'F') || (T[i] == 'V') ) && ((T[i-1] != 'B') && (T[i-1] != 'P') && (T[i-1] != 'F') && (T[i-1] != 'V'))) S += '1';
   else if (((T[i] == 'C') || (T[i] == 'S') || (T[i] == 'K') || (T[i] == 'G') || (T[i] == 'J') || (T[i] == 'Q') || (T[i] == 'X') || (T[i] == 'Z') ) && ((T[i-1] != 'C') && (T[i-1] != 'S') && (T[i-1] != 'K') && (T[i-1] != 'G') && (T[i-1] != 'J') && (T[i-1] != 'Q') && (T[i-1] != 'X')
               && (T[i-1] != 'Z') ) ) S += '2';
   else if (((T[i] == 'D') || (T[i] == 'T') ) && ((T[i-1] != 'D') && (T[i-1] != 'T') ) ) S += '3';
   else if ((T[i] == 'L') && (T[i-1] != 'L') ) S += '4';
   else if (((T[i] == 'M') || (T[i] == 'N') ) && ((T[i-1] != 'M') && (T[i-1] != 'N'))) S += '5';
   else if ((T[i] == 'R') && (T[i-1] != 'R')) S += '6';
        m = S.size();
    }
    for(int i = 4; i > m; i--)
    {
        S += '0';
    }
    return S;
}

int main(){

   cout<<"         NAME                     SOUNDEX CODE"<<endl;

while(cin >> T)
{

Soundex = geraSoundex(T);
while(T.size() != 25)
  T +=  ' ';
cout <<"         " << T << Soundex << endl;

}

cout << "                   END OF OUTPUT" << endl;
}
