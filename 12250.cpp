//Anderson Zudio
//AC em 3/7/2014

#include <iostream>

using namespace std;

int main(){
    int t = 0;
    string x;
    while(cin >> x, x != "#"){
        if(x == "HELLO") cout << "Case " << ++t << ": ENGLISH" << endl;
        else if(x == "HOLA") cout << "Case " << ++t << ": SPANISH" << endl;
        else if(x == "HALLO") cout << "Case " << ++t << ": GERMAN" << endl;
        else if(x == "BONJOUR") cout << "Case " << ++t << ": FRENCH" << endl;
        else if(x == "CIAO") cout << "Case " << ++t << ": ITALIAN" << endl;
        else if(x == "ZDRAVSTVUJTE") cout << "Case " << ++t << ": RUSSIAN" << endl;
        else cout << "Case " << ++t << ": UNKNOWN" << endl;
    }
    return 0;
}
