//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 05/04/2014
#include <iostream>
using namespace std;
int a, b, n;

int PotMod(int a, int b, int n){
    long long int m, c, d, rb[101];  int i, j;
    m= b;   i= 101;
    while (m > 0){
        rb[--i]= m % 2;   m= m /2;
    }
    c= 0; d= 1;
    for(j=i; j<= 100; j++){
        d= (d*d) % n;    c= 2*c;
        if (rb[j] == 1){
            d= (d*a) % n;   c++;
        }
        //cout << rb[j] << " " << c << " " << d <<endl;
    }
    return (int) d;
}

int main(){
    while(cin >> a >> b >> n){

        cout << PotMod(a, b, n)<<endl;
    }
    return 0;
}
