//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 29/03/2014

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <math.h>

using namespace std;

typedef struct ret { long long int d; long long int x; long long int y;};
int a, b; ret r;

ret MDCE(long long int a, long long int b, ret ra)
{   long long int xx;
    if (b == 0)
    {   ra.d = a;  ra.x = 1;  ra.y = 0;
    }
    else
    {   ra = MDCE(b, a % b, ra);
        xx = ra.x;  ra.x = ra.y;  ra.y = xx - (a / b)*ra.y;
    }
    return ra;
}

int main()
{
    long long int a, b;
    while(cin >> a >> b)
    {
        r = MDCE(a, b, r);
        cout << r.x << " " << r.y << " " << r.d << endl;
    }

}
