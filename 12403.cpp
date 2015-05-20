//Anderson Zudio
//AC em 1/7/2014

#include<iostream>

using namespace std;

int main()
{
    int t, resp, y;
    string x;
    cin >> t;
    resp = 0;
    for(int i = 0; i < t; i++){
        cin >> x;
        if(x == "report") cout << resp << endl;
        else{
            cin >> y;
            resp += y;
        }
    }
    return 0;

}
