//Anderson Zudio
//AC em 3/7/2014

#include <iostream>

 using namespace std;

int main()
{
    int t, h, m, time;
    char lixo;
    cin >> t;
    while (t--)
    {
        cin >> h >> lixo >> m;;
        time = h*60 + m;
        if (time > 720) time -= 720;
        time = 720 - time;

        h = (time/60);
        m = time%60;

        if (!h) h = 12;
        h < 10 ? cout << '0' << h << ':' : cout << h << ':';
        m < 10 ? cout << '0' << m << endl : cout << m << endl;
    }
    return 0;
}
