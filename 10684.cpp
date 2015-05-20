//AC em 26/04/204
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar

#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;
long long int S[100050], n;





void SCSM(){
    long long int sm, ss;
 sm=0;
 ss=0;
 for (int i=1; i<=n; i++){
  if ((ss + S[i]) >= 0){
   ss=ss + S[i];
   if (ss > sm)
                sm=ss;

        }
  else
   ss=0;

 }

 if (sm>0)
 cout <<"The maximum winning streak is "<< sm<<".";
 else
 cout<< "Losing streak.";

}




int main()
{

while(cin>>n, n != 0){
    for (int i=1; i<=n; i++) cin >> S[i];
    SCSM();
    cout << endl;
}



return 0;
}
