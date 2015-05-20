//Anderson Zudio
//Victor Cracel
//Julio Cesar
//Ac em 3/9/2014
#include <iostream>
#include <queue>
#include <vector>  
#include <functional> 
using namespace std;

vector <pair <int, int> > v;
vector <int> v2;
priority_queue<int> p;
int n, j;
pair <int, int> x;

int main() {
	while(cin >> n, n){
		v.clear();
		v2.clear();
		for(int i = 0; i < n; i++){
			cin >> x.first >> x.second;
			v.push_back(x);
		}
		p.push(v[n-1].first);
		v2.push_back(v[0].first);
		auto bola = v[0].second;
		for(int i = n-2; i > 0; i--){
			bola += v[i].second;
			p.push(v[i].first);
		}
		
		while(bola && !p.empty()){ 
			v2.push_back(p.top()/2);
			p.pop();
			bola--;
		}
		
		auto resp = 0;
		for(auto i = v2.begin(); i != v2.end(); i++){
			resp += *i;
		}
		
		while(!p.empty()){
			resp += p.top();
			p.pop();
		}
		cout << resp << endl;
		
	}
	return 0;
}