//AC em 24/05/2014
//Anderson Zudio de Moraes 
//Julio Cesar
//Victor Cracel Messner
#include <iostream> 
#include <cmath> 
#define PI acos(-1) 
  
double chord, arc, s, a; 
char ang[3]; 
  
using namespace std; 
  
int main() 
{ 
    cout.setf(ios::fixed); 
    cout.precision(6); 
  
    while(cin >> s >> a >> ang) 
    { 
    if (ang[0] == 'd') a *= PI/180; 
    else a *= PI/10800; 
    s += 6440; 
    a = acos(cos(a)); 
    arc = a*s; 
    chord = 2*sin(a/2)*s; 
  
    cout << arc << " " << chord << endl; 
    } 
  
}