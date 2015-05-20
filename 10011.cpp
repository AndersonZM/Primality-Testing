//Anderson Zudio de Moraes
//Victor Cracel Messner
//Julio Cesar
//AC em 30/05/2014
#include<iostream>
#include<cmath>

#define PI acos(-1)

using namespace std;

const long double QZERO = 0.000000000000001;

typedef struct { long double x, y; }ponto;
typedef struct { long double a,b,c; }reta;
typedef struct {ponto pc; long double raio; }circulo;
typedef struct { reta r1,r2; }duasretas;
typedef struct { ponto p1,p2; }doispontos;


double seguro_sqrt(double x) {
    if (fabs(x) < QZERO) {
        return 0;
    }
    else{
        return sqrt(x);
    }
}


long double DistPontos(ponto p1,ponto p2){
    return (seguro_sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)));
};

long double AnguloPontoOrigem(ponto p1){
    long double ang;
    if (fabs(p1.x)<QZERO){
        if (fabs(p1.y)< QZERO) ang = 0;
        else if (p1.y < 0)     ang = 3*PI/2;
        else                   ang = PI/2;
    }
    else{
        ang = atan(p1.y/p1.x);
        if (p1.x < 0)       ang = ang + PI;
        else if (p1.y < 0)  ang = ang + 2*PI;
    }
    return ang;
}

ponto PontoIntersecao(reta r1,reta r2){
    ponto p;
        p.x = (r2.b*r1.c-r1.b*r2.c)/(r2.a*r1.b-r1.a*r2.b);
        if (fabs(r1.b)>QZERO)
            p.y = -(r1.a*p.x + r1.c)/r1.b;
        else
            p.y = -(r2.a*p.x + r2.c)/r2.b;
    return p;
};

reta PontosParaReta(ponto p1,ponto p2){
    reta r;
    if (p1.x==p2.x){ r.a=1; r.b =0; r.c=-p1.x;}
    else{ r.b=1; r.a=-(p1.y-p2.y)/(p1.x-p2.x);
          r.c=-(r.a*p1.x)-(r.b*p1.y);
    }
    return r;
};

reta PontoTangenteParaReta(ponto p, long double tang){
    reta r;
    r.a=-tang;  r.b=1;  r.c=-(r.a*p.x + r.b*p.y);
    return r;
};

reta RetaParalela(ponto p1, reta r1){
    reta r2;
    r2.a = r1.a;  r2.b = r1.b;  r2.c = -(r1.a*p1.x + r1.b*p1.y);
    return r2;
};

reta RetaPerpendicular(ponto p1, reta r1){
    ponto p2; reta r2;
    if (fabs(r1.b)<= QZERO){
        r2.a = 0;  r2.b = 1;  r2.c = -p1.y;
    }
    else if (fabs(r1.a)<= QZERO){
        r2.a = 1;  r2.b = 0;  r2.c = -p1.x;
    }
    else r2 = PontoTangenteParaReta(p1,1/r1.a);
    return r2;
};

ponto PontoMaisProximo(ponto p1, reta r1){
    ponto p2;
    reta r2;
    if (fabs(r1.b)<= QZERO){
        p2.x = -r1.c;  p2.y = p1.y;
    }
    else if (fabs(r1.a)<= QZERO){
        p2.x = p1.x;  p2.y = -r1.c;
    }
    else{
        r2 = PontoTangenteParaReta(p1,1/r1.a);
        p2 = PontoIntersecao(r1,r2);
    }
    return p2;
};

ponto RotacaoHoraria(ponto p, long double ang/*em radianos*/){
    ponto p1;
    p1.x = p.x*cos(ang)-p.y*sin(ang);
    p1.y = p.x*sin(ang)+p.y*cos(ang);
    return p1;
};

duasretas Tangentes(ponto p1, circulo c){
    duasretas dr; long double alfa,beta,d,b; ponto p0; circulo c1; int nrot;
    c1.raio = c.raio;  c1.pc.x = c.pc.x-p1.x; c1.pc.y = c.pc.y-p1.y;
    p0.x =0;  p0.y = 0;  nrot = 0;
	while (((c1.pc.x < 0) && (fabs(c1.pc.x) > QZERO)) ||
          ((c1.pc.y < 0) && (fabs(c1.pc.y) > QZERO))) {
        nrot++;  c1.pc = RotacaoHoraria(c1.pc, PI/2);
    }
 	d = DistPontos(p0,c1.pc);   b = seguro_sqrt(d*d-c1.raio*c1.raio);
    alfa = atan(c1.raio/b);
    if (nrot > 1) nrot = nrot -2; // Mesma tangente
	if (fabs(c1.pc.x) <= QZERO) beta = PI/2*(nrot-1);
	else                        beta = atan(c1.pc.y/c1.pc.x) -PI/2*nrot;
    if (fabs(cos(beta-alfa)) <= QZERO) {
        dr.r1.a = 1;  dr.r1.b = 0;  dr.r1.c = -p1.x;
    }
    else dr.r1 = PontoTangenteParaReta(p1,sin(beta-alfa)/cos(beta-alfa));
    if (fabs(cos(beta+alfa)) <= QZERO) {
        dr.r2.a = 1;  dr.r2.b = 0;  dr.r2.c = -p1.x;
    }
    else dr.r2 = PontoTangenteParaReta(p1,sin(beta+alfa)/cos(beta+alfa));
	return dr;
}

doispontos PontosTangencia(ponto p1, circulo c){
    duasretas dr;  doispontos dp;
    dr = Tangentes(p1,c);
    dp.p1 = PontoMaisProximo(c.pc,dr.r1);
    dp.p2 = PontoMaisProximo(c.pc,dr.r2);
    return dp;
}

void Solucao(long double y, long double x, long double r, long double u, long double v){
    ponto p0, ph, p1, p2; circulo c;  duasretas dr; doispontos dp; reta tangente1, tangente2;
    long double angh, angt1, angt2, resp;

    c.pc.x = x; c.pc.y = y; c.raio = r; ph.x = v; ph.y = u;
    p0.x = p0.y = 0;

    dr = Tangentes(p0,c);
    tangente1 = dr.r1;
    tangente2 = dr.r2;

    dp = PontosTangencia(p0,c);
    p1 = dp.p1;
    p2 = dp.p2;

    //cout<<"Primeira tangente: "<<dr.r1.a<<"x + "<<dr.r1.b<<"y + "<<dr.r1.c<<" = 0"<<endl;
    //cout<<"Segunda  tangente: "<<dr.r2.a<<"x + "<<dr.r2.b<<"y + "<<dr.r2.c<<" = 0"<<endl;

    angh = AnguloPontoOrigem(ph); angh = 180*angh / PI;
    angt1 = AnguloPontoOrigem(p1); angt1 = 180*angt1 / PI;
    angt2 = AnguloPontoOrigem(p2); angt2 = 180*angt2 / PI;
    if(angt1 == angt2){
        c.raio -= 0.0000000001;
        dp = PontosTangencia(p0,c);
        p1 = dp.p1;
        p2 = dp.p2;
        dr = Tangentes(p0,c);
        tangente1 = dr.r1;
        tangente2 = dr.r2;
        angh = AnguloPontoOrigem(ph); angh = 180*angh / PI;
        angt1 = AnguloPontoOrigem(p1); angt1 = 180*angt1 / PI;
        angt2 = AnguloPontoOrigem(p2); angt2 = 180*angt2 / PI;
        c.raio += 0.0000000001;
    }
     //      cout << p1.x << p1.y << endl;
    resp = 0;
  // cout << angh << " " << angt1 << " " << angt2 << endl;
     if((angt1 <= angh) && (angt2 >= angh)){
           // cout << DistPontos(c.pc,p0) << " " << DistPontos(ph,p0) << endl;
            if((DistPontos(p1,p0) <= DistPontos(ph,p0)) && (DistPontos(p2,p0) <= DistPontos(ph,p0))){
                p1 = PontoMaisProximo(ph,tangente1);
                p2 = PontoMaisProximo(ph,tangente2);
                //cout << DistPontos(ph,p2) << " " << DistPontos(ph,p1) << " " << DistPontos(ph,c.pc) - c.raio << endl;
                resp = min(DistPontos(ph,p2),DistPontos(ph,p1));
                resp = min(resp,DistPontos(ph,c.pc) - c.raio);

            }
    }
    cout << fabs(resp) << endl;
}

int main(){
    int t;
    long double x, y, r, u, v;
    cin >> t;
    cout.setf(ios_base::fixed);
    cout.precision(3);
    while(t--){
            cin >> y >> x >> r >> u >> v;
            if(x > 0){
                x = -x;
                v = -v;
            }
            Solucao(y, x, r, u, v);
    }
    return 0;
}

