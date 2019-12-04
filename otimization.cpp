#include <iostream>
#include <vector>

using namespace std;
using matrix = vector<vector<float>>;

// Função

float xb[] = {x1,x2};

float[] mult(*x1, *x2)
{
	int x1_cols = sizeof(x1[0])/sizeof(x1[0][0]);
  int x2_rows = sizeof(x2)/sizeof(x2[0]);
	if (x1_cols == x1_rows){
			int c[a_rows][b_cols];
			int i, j,k;
			for (i=0;i<a_rows; i++){
					for (j=0;j<b_cols;j++){
							c[i][j]=0;
							for (k=0;k<b_rows;k++){
									c[i][j] = c[i][j]+a[i][k]*b[k][j];
							}
					}
			}
	}else{
		return {{}}
	}
}

void func(x1,x2) {
	float res = -log((x2*x1)/(1+(x1*x2)^2));
	return res;
}

void phi(x1,x2,dx,dy,t) {
	return func(x1+t*dx, x2+t*dy);
}

void gradf(x1,x2) {
	vector<float> res;
	res.push_back(((x1^2)*(x2^2)-1)/((x1^3)*(x2^2)+x1));
	res.push_back(((x1^2)*(x2^2)-1)/((x1^2)*(x2^3)+x2));
	return res;
}

// Seção Áurea
void GoldenSection(x1,x2){
	float eps = 1, ro = 1;
	float dt1 = (3-sqrt(5))/2;
	float dt2 = 1-dt1;
	float a = 0, s = ro , b = 2*ro;
	
	
	while (phi(x1,x2,d1,d2,b) < phi(x1,x2,d1,d2,s){
		a = s;
		s = b;
		b = 2*b;
	}
	
	float 	u = a+dt1*(b-a), 
			v = a+dt2*(b-a),
			t = 1;
	
	while ((b-a)>eps){
		if (phi(x1,x2,d1,d2,u)<phi(x1,x2,d1,d2,v)){
			b=v, v=u, u= (a+dt1(b-a)); // também tem que arrumar, é array
		}
		else {
			a = u, u = v, v = (a+dt2(b-a));
		}
		return t = (u+v)/2; // isso é inteiro??? mas tem array .-. 
	}
}

// Armijo

void Armijo(x1,x2){
	vector<float> d;
	float y, n, t=1;
	float grad[] = gradf(x1,x2);
	while ( phi(x1,x2,d.at(0),d.at(1),t) > func(x1,x2)+ n*t*(grad[0]*d.at(0) + grad[1]*d.at(1)))
	{
		t = y*t;
	}
}

// Método do Gradiente

void Gradient(x1,x2)[
	vector<float> grad = gradf(x1,x2);
	vector<float> x = x.push_back(x1).push_back(x2);
	vector<float> xk;
	if (grad.at(0) != 0 && grad.at(1) != 0) {
		vector<float> d = d.push_back(-grad.at(0)).push_back(-grad.at(1));
		// obter t>0 onde f(x+td)>f(x) por Armijo ou Seção Áurea
		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
		return Gradient(xk.at(0),xk.at(1));
	}
	else {
		return x;
	}
]

// Método de Newton

void Newton(x1,x2) {
	vector<float> grad = gradf(x1,x2);
	vector<float> x = x.push_back(x1).push_back(x2);
	vector<float> xk;
	if (grad.at(0) != 0 && grad.at(1) != 0) {
		vector<float> d = // calcular a hessiana inversa
		//obter t>0
		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
		return Newton(xk.at(0),xk.at(1));
	}
	else {
		return x[]
	}
}

// Método Quase-Newton

void QuasiNewlton(x1,x2,hess) {
	vector<float> grad = gradf(x1,x2);
	vector<float> x = x.push_back(x1).push_back(x2);
	vector<float> xk;
	matrix H = hess;
	if (grad.at(0) != 0 && grad.at(1) != 0) {
		vector<float> d = d.push_back(-(H.at(0)*grad.at(0)+H.at(1)*grad.at(1)).push_back(-(H.at(2)*grad.at(0)+H.at(3)*grad.at(1))
		//obter t>0
		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
		// determinar proximo H
		vector<float> gradk = gradf(xk.at(0),xk.at(1));
		vector<float> p = p.push_back(xk.at(0) - x.at(0)).push_back(xk.at(1) - x.at(1));
		vector<float> q = q.push_back(gradk.at(0) - grad.at(0).push_back(gradk.at(1) - grad.at(1));
		matrix Hk;
		return QuasiNewlton(xk.at(0),xk.at(1));
	}	else {
		return x[];
	}
}