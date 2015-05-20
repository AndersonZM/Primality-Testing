//AC 24/05/2014
//Anderson Zudio de Moraes 
//Victor Cracel Messner
//Julio Cesar
#include<iostream> 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#define PI acos(-1) 
  
double a,b,c, s, A1,r,R,A2,A3,Sun,Vio,Ros; 
  
using namespace std; 
  
int main() 
{ 
	while(cin >> a >> b >> c){

    	s = a+b+c; 
		s /= 2; 
	    A1 = sqrt(s*(s-a)*(s-b)*(s-c)); 
	    r = A1/s; 
	    R = (a*b*c)/(4*A1); 
		A2 = PI*r*r; 
    	A3 = PI*R*R; 
	    Sun = A3-A1; 
		Vio = A1-A2; 
    	Ros = A2; 
    	cout.setf(ios::fixed); 
   		cout.precision(4); 
    	cout << Sun << " " << Vio << " " << Ros << endl; 
	}
  
} 
