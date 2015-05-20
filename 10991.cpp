//AC em 24/05/2014
//Anderson Zudio de Moraes 
//Victor Cracel Messner
//Julio Cesar
 
#include <iostream>
#include <math.h>
#define PI acos(-1)
using namespace std;
 
double r1, r2, r3, l1, l2 , l3, atri, x, y, z, s1, s2, s3, res;
int t;
 
int main() {
 
	cin >> t;
	while(t--){
		cin >> r1 >> r2 >> r3;
 
		l1 = r1 + r2;
		l2 = r1 + r3;
		l3 = r2 + r3;
 
		//cout << l1 << " " << l2 << " " << l3 << endl;
 
		x = acos(  (l2*l2 - l1*l1 - l3*l3)  /(((-2) * l1 * l3)));
		y = acos(  (l1*l1 - l2*l2 - l3*l3)  /(((-2) * l2 * l3)));
		z = acos(  (l3*l3 - l1*l1 - l2*l2)  /(((-2) * l1 * l2)));
 
		//cout << x << " " << y << " " << z << endl;
 
		s1 = (r1*r1*z) / 2.0;
		s2 = (r2*r2*x) / 2.0;
		s3 = (r3*r3*y) / 2.0;
 
		//cout << s1 << " " << s2 << " " << s3 << endl;
 
 
		atri = (l1 * l2 * sin(z)) / 2.0;
 
	cout.setf(ios::fixed); 
     	cout.precision(6);
     	res = atri - s1 - s2 - s3;
		cout << res << endl;
 
	}
 
	return 0;
}