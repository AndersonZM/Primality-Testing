//AC em 24/05/2014
//Anderson Zudio de Moraes 
//Victor Cracel Messner
//Julio Cesar
#include <iostream>
#include <cmath>
#define PI acos(-1) 
  
using namespace std; 
  
char cidov; //Coisa inutil  do valladolid 
double comp,larg, r, ang; 
  
int main() 
{ 
   int tcase = 0; 
    cout.setf(ios_base::fixed);
    cout.precision(10); 
    while(cin >> comp >> cidov >> larg) 
    { 
       r = sqrt((comp*comp/4)+(larg*larg/4)); 
       ang = asin((larg/2)/r); 
       larg *= 400/((2*comp)+(4*ang*r)); 
       comp *= 400/((2*comp)+(4*ang*r)); 
  
       cout << "Case " << ++tcase << ": " << comp <<  " " << larg << endl; 
    } 
} 