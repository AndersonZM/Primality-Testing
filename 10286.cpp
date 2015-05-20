//AC em 01/02/2013
//Anderson Zudio de Moraes

#include <iostream>
#include <cmath>
#define PI acos(-1)

using namespace std;

double F, h, a, b;

int main()
{

    cout.setf(ios_base::fixed);
    cout.precision(10);

    while(cin >> F)
    {
       h = sin(9*PI/180) * F;
       cout << (h/tan(9*PI/180)) + (h/tan(63*PI/180)) << endl;
    }
}
