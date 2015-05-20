//Anderson Zudio
//AC em 30/07/2014
#include <iostream>
#include<queue>
#include <vector>
using namespace std;

int main() {
	int n, x, resp, s;
	auto comp = [](int a, int b)->bool{ return a > b; };
	priority_queue <int, vector<int>, decltype(comp) > queue(comp);
	while(cin >> n, n){
		for(int i = 0; i < n; i++){
			cin >> x;
			queue.push(x);
		}
		if(n == 1){ cout << x << endl; continue; }
		s = 0;
		while(queue.size() > 1){
            x = queue.top();
            queue.pop();
            x += queue.top();
            queue.pop();
            s += x;
            queue.push(x);
		}
		queue.pop();
		cout << s << endl;
	}
	return 0;
}
