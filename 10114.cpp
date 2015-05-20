//Anderson Zudio
//AC em 31/07/2014
#include<iostream>
#include<cmath>
#include <vector>

using namespace std;


int main()
{
	int d, n, m;
	double a, loan, aux;
	vector <double> hist(1010);
	while(cin >> d){
		if(d < 0) break;
		double a, loan;
		cin >> a >> loan >> n;
		fill(hist.begin(), hist.begin()+1000, 0);
		for(int i = 0; i < n;i++){
			cin >> m >> aux;
			hist[m] = aux;
		}
		auto c = (double)(a + loan) * (1.0 - hist[0]);
		double loan_v = loan / (double) d;
		if(c > loan){
			cout << "0 months" << endl;
			continue;
		}
		for(int i = 1 ; i <= d ; i ++)
		{
			loan -= loan_v;
			double *de_value = &hist[i];
			while(fabs(*de_value-0) < 1e-6) de_value--;
			c *= (1.0 - *de_value);


			if(c > loan){
				if(i > 1)
					cout << i << " months" << endl;
				else
					cout << "1 month" << endl;
				break;
			}
		}
	}
	return 0;
}
