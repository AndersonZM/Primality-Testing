//AC 26/04/2014
//Anderson Zudio
//Julio Cesar

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

int i,p,n,k,t;
string s;

void resolveAi() {
    bool b = true;
    int metade=pow(2,n)/2;
    while(metade) {
        if (b && k <= metade) { s += "0"; }
        else if (b && k > metade) { s += "1"; k-=metade; b=false; }
        else if (!b && k <= metade) { s+= "1"; b=true; }
        else if (!b && k > metade) { s += "0"; k-=metade; }
        metade/=2;
    }
}

int main() {
    long long int resp;
    cin>>t;
    while(t) {
        cin>>n>>k;
        ++k;
        s="";
        resolveAi();
        resp=0; p=0;
        for(i=s.length()-1;i>=0;--i,++p) {
            if (s[i]=='1')
                resp+=pow(2,p);
        }
        cout<<resp<<endl;
        --t;
    }
    return 0;
}
