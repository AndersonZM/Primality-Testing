//Acha os dois pontos com menor distancia entre si em um conjunto de pontos
//C++11 implementação O(n log(n))
//Anderson Zudio
//Victor Cracel Messner
//Julio Cesar
//AC 09/06/2014

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#define QZERO 0.000001

using namespace std;

typedef struct ponto{
    double x;
    double y;
}Ponto;

typedef struct doispontos{
    Ponto p1;
    Ponto p2;
}DoisPontos;

double DistanciaEntrePontos(Ponto a, Ponto b){
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

//Solução força bruta O(n^2), pode ser usada para qualquer quantidade de pontos
pair<double, DoisPontos> forcBruta(vector<Ponto> Pontos){
	if (Pontos.size() < 2){
        DoisPontos ParPonto;
        ParPonto.p1.x = ParPonto.p1.y = ParPonto.p2.x = ParPonto.p2.y = 0;
        auto resp = pair<double, DoisPontos>(-1, ParPonto);
		return resp;
	}
	auto menorDist = abs(DistanciaEntrePontos(Pontos.at(0), Pontos.at(1)));
	DoisPontos menorPar;
	menorPar.p1 = Pontos.at(0); menorPar.p2 = Pontos.at(1);
	for (auto i = begin(Pontos); i != (end(Pontos) - 1); i++){
		for (auto j = i+1; j < end(Pontos); j++){
			auto dist = abs(DistanciaEntrePontos(*i, *j));
			if (dist < menorDist){
				menorDist = dist;
				menorPar.p1 = *i;
				menorPar.p2 = *j;
			}
		}
	}
    auto resp = pair<double, DoisPontos>(menorDist, menorPar);
    return resp;
}

//Solução O(n log (n)) para n > 3
pair<double, DoisPontos> MenorDistParPonto(vector<Ponto> PontosX, vector<Ponto> PontosY){
    auto n = PontosX.size();
	if (n <= 3) {
		return forcBruta(PontosX);
	}
	auto PontosXl = vector<Ponto>();
	auto PontosXr = vector<Ponto>();
	copy(begin(PontosX), begin(PontosX) + (n/2), back_inserter(PontosXl));
	copy(begin(PontosX) + (n/2), end(PontosX), back_inserter(PontosXr));
	auto xMed = PontosX.at(n/2).x;
	auto PontosYl = vector<Ponto>();
	auto PontosYr = vector<Ponto>();
	copy_if(begin(PontosY), end(PontosY), back_inserter(PontosYl), [&xMed](Ponto p) { return p.x <= xMed; });
	copy_if(begin(PontosY), end(PontosY), back_inserter(PontosYr), [&xMed](Ponto p) { return p.x > xMed; });
	auto L = MenorDistParPonto(PontosXl, PontosYl);
	auto R = MenorDistParPonto(PontosXr, PontosYr);
	auto menorPar = (L.first <= R.first) ? L : R;
	auto PontosYs = vector<Ponto>();
	copy_if(begin(PontosY), end(PontosY), back_inserter(PontosYs), [&menorPar, &xMed](Ponto p) { return abs(xMed - p.x) < menorPar.first;	});
	auto resp = menorPar;
	auto tYs = PontosYs.size();
	for (auto i = 0; i < tYs-1; i++) {
		for (auto j = i+1; j < tYs && ((PontosYs.at(j).y - PontosYs.at(i).y) < menorPar.first); j++) {
			auto dist = abs(DistanciaEntrePontos(PontosYs.at(j), PontosYs.at(i)));
			if (dist < resp.first) {
                DoisPontos ParPonto;
                ParPonto.p1 = PontosYs.at(j);
                ParPonto.p2 = PontosYs.at(i);
				resp = { dist, ParPonto };
			}
		}
	}
	return resp;
}

//Impressão da resposta
void imprimePontos(pair<double, DoisPontos> resposta){
        double distancia;
        Ponto p1, p2;
        //Distancia guarda a distancia do par de pontos da resposta
        distancia = resposta.first;
        //Variavel p1 e p2 recebem os pontos de resposta
        p1 = resposta.second.p1;
        p2 = resposta.second.p2;
        if(distancia >= 0 && distancia < 10000)
            cout << distancia << endl;
        else
            cout << "INFINITY" << endl;
}

//Busca solução trivial (dois pontos iguais)
bool verifica(vector<Ponto> PontosX){
    auto achei = false;
        for(auto i = begin(PontosX), j = i+1; i != end(PontosX) - 1; i++, j++)
            if(i->x == j->x && i->y == j-> y){
                DoisPontos pontoIgual;
                pontoIgual.p2 = pontoIgual.p1 = *i;
                pair<double, DoisPontos> trivial = { 0, pontoIgual};
                imprimePontos(trivial);
                achei = true;
                break;
            }
            return achei;
}

int main() {
    //Declaracões gerais
    vector<Ponto> Pontos;
    int np;
    cout.setf(ios::fixed);
    cout.precision(4);

    while(cin >> np, np){

        Pontos.clear();
        for(int i = 0; i < np; i++){
            Ponto p;
            cin >> p.x >> p.y;
            Pontos.push_back(p);
        }

        //Ordenação de pontos por X
        sort(begin(Pontos), end(Pontos), [](Ponto a, Ponto b) { if(a.x != b.x) return a.x < b.x; else return a.y < b.y; });
        auto PontosX = Pontos;
        if(verifica(PontosX)) continue;

        //Ordenação de pontos por Y
        sort(begin(Pontos), end(Pontos), [](Ponto a, Ponto b) { return a.y < b.y; });
        auto PontosY = Pontos;

        //Resposta O(n^2) força bruta se necessário
        /*auto resposta = forcBruta(PontosX);*/

        //Resposta em O(nlog(n))
        auto resposta = MenorDistParPonto(PontosX, PontosY);
        imprimePontos(resposta);
    }
	return 0;
}
