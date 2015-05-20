//Anderson Zudio
//AC em 3/7/2014

#include<iostream>

using namespace std;

int main()
{
    int t = 0;
    string x;
    while(cin >> x, x != "*"){
        if(x == "Hajj") cout << "Case " << ++t << ": Hajj-e-Akbar" << endl;
        else cout << "Case " << ++t << ": Hajj-e-Asghar" << endl;
    }
    return 0;

}
