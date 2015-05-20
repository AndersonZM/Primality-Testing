//Anderson Zudio
//Ac em 2/7/2014

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(){
    string line, word;
    bool what = true;
    while(getline(cin, line)){
            for(auto i = 0; i < line.size(); i++){
                if(line[i] == '"'){
                    what ? line.insert(i, "``") : line.insert(i, "''");
                    what = !what;
                    line.erase(i+2, 1);
                }
            }
            cout << line <<endl;
        }
    return 0;
}
