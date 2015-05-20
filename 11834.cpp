//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 24/05/2014

#include <iostream>
using namespace std;
int main(){
    while (true){
        int l, c, r1, r2,hip1,hip2,hip3,hip4, diagonal,total,maiorr,menorl;
        cin >> l >> c >> r1 >> r2;
        if (!(l && c && r1 && r2)) {break;} 
         
        if (l<c){menorl=l;}
        else {menorl=c;}
        if (r1>r2){maiorr=r1*2;}
        else{maiorr=r2*2;}
         
        hip4= (l-(r1+r2))*(l-(r1+r2)) + (c-(r1+r2))*(c-(r1+r2)) ;
        if((hip4>=(r1+r2)*(r1+r2)) && (maiorr <= menorl)){cout<<"S"<<endl;}
        else {cout << "N" << endl;}}
    return 0;}