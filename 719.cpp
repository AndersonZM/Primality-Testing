//Ac em 10/05/2014 
//Anderson Zudio 
//Victor Cracel Messner 
//Julio Cesar  
#include <stdio.h>
#include <string.h>
#include<iostream>
#define MAX 11000
 
using namespace std;
 
int pior(char *str, int len) {
    int i = 0, j = 1, k = 0;
    int a, b;
    while(i < len && j < len && k < len) {
        a = (i+k >= len) ? i+k-len : i+k;
        b = (j+k >= len) ? j+k-len : j+k;
        if(str[a] == str[b])
            k++;
        else {
            str[a] > str[b] ? i += k+1 : j += k+1;
            if(i == j)
                j++;
            k = 0;
        }
    }
 
    return i;
}
 
int main() {
    int T;
    char str[MAX];
    cin >> T;
    while(T--) {
        cin >> str;
        cout << pior(str, strlen(str))+1 << endl;
    }
    return 0;
}