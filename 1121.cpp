	//Anderson Zudio
	//AC 02/09/2014
    #include <iostream>
    #include <vector>
    #define INF 1000000
    using namespace std;
     
    vector <int> vetor(300000);
     
    int main() {
     
    int n, a;
     
    while(cin >> n >> a){
    
		for(int i = 0; i < n; i++)
		  	cin >> vetor[i];
    	
    	auto soma = vetor[0];
    	auto resp = INF;
    	
	 	for(int i = 0, j = 0; i < n;) {
			 if (soma >= a) {
	    		auto aresp = i - j + 1;
	    		resp = min(aresp, resp);
	    		if(j < i) soma -= vetor[j++];
	    	}
	    	else{
	    		i++;
	    		if(i < n) soma += vetor[i];
	    	}
    	}
    	cout << (resp==INF ? 0:resp) << endl;
    }
  
    return 0;
}