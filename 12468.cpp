//Anderson Zudio
//AC em 31/07/2014
#include <iostream>

using namespace std;

int main(){
    int ch1, ch2;
    while(cin >> ch1 >> ch2, ch1 != -1){
        int aux = ch1, first = 0, second = 0;
        while(aux != ch2){
            first++;
            aux++;
            aux = aux % 100;
        }
        aux = ch1;
        while(aux != ch2){
            second++;
            if(aux == 0) aux = 99;
            else aux--;
        }
        cout << min(first, second) << endl;
    }
    return 0;
}
