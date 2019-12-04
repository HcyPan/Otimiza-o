#include <iostream>
#include <vector>

using namespace std;
using matrix = vector<vector<float> > ;

// Função
class MathFunctions {
  public: 
    enum getStepFunction { Armijo = 1, GoldenSect = 2 };

    void printMat(const matrix &C);
    void printVec(const vector<float> &C);
    
    void initializeMat(matrix &A);
    void initializeVec(vector<float> &A);

    matrix matMul(const matrix &A, const matrix &B);
    vector<float> matVecMul(const vector<float> &A, const matrix &B);
    float vecMul(const vector<float> &A, const vector<float> &B);
    void vecMulFloat(vector<float> &A, float B);
    matrix matSum(const matrix &A, const matrix &B);
    matrix matSumValue(const matrix &A, const float B);
    
    float funcObj(const vector<float> &X);
    float phi(const vector<float> &X, const vector<float> &dir, float t);
    vector<float> gradF(const vector<float> &X);
    float goldenSection(const vector<float> &X,  const vector<float> &dir);
    float armijo(const vector<float> &X,  const vector<float> &dir);
    vector<float> gradient(const vector<float> &X, int stepFunction);
};

// // Seção Áurea
// void GoldenSection(x1,x2){
// 	float eps = 1, ro = 1;
// 	float dt1 = (3-sqrt(5))/2;
// 	float dt2 = 1-dt1;
// 	float a = 0, s = ro , b = 2*ro;bele
	
	
// 	while (phi(x1,x2,d1,d2,b) < phi(x1,x2,d1,d2,s){
// 		a = s;
// 		s = b;
// 		b = 2*b;
// 	}
	
// 	float 	u = a+dt1*(b-a), 
// 			v = a+dt2*(b-a),
// 			t = 1;
	
// 	while ((b-a)>eps){
// 		if (phi(x1,x2,d1,d2,u)<phi(x1,x2,d1,d2,v)){
// 			b=v, v=u, u= (a+dt1(b-a)); // também tem que arrumar, é array
// 		}
// 		else {
// 			a = u, u = v, v = (a+dt2(b-a));
// 		}
// 		return t = (u+v)/2; // isso é inteiro??? mas tem array .-. 
// 	}
// }

// // Armijo

// void Armijo(x1,x2){
// 	vector<float> d;
// 	float y, n, t=1;
// 	float grad[] = gradf(x1,x2);
// 	while ( phi(x1,x2,d.at(0),d.at(1),t) > func(x1,x2)+ n*t*(grad[0]*d.at(0) + grad[1]*d.at(1)))
// 	{
// 		t = y*t;
// 	}
// }

// // Método do Gradiente

// void Gradient(x1,x2)[
// 	vector<float> grad = gradf(x1,x2);
// 	vector<float> x = x.push_back(x1).push_back(x2);
// 	vector<float> xk;
// 	if (grad.at(0) != 0 && grad.at(1) != 0) {
// 		vector<float> d = d.push_back(-grad.at(0)).push_back(-grad.at(1));
// 		// obter t>0 onde f(x+td)>f(x) por Armijo ou Seção Áurea
// 		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
// 		return Gradient(xk.at(0),xk.at(1));
// 	}
// 	else {
// 		return x;
// 	}
// ]

// // Método de Newton

// void Newton(x1,x2) {
// 	vector<float> grad = gradf(x1,x2);
// 	vector<float> x = x.push_back(x1).push_back(x2);
// 	vector<float> xk;
// 	if (grad.at(0) != 0 && grad.at(1) != 0) {
// 		vector<float> d = // calcular a hessiana inversa
// 		//obter t>0
// 		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
// 		return Newton(xk.at(0),xk.at(1));
// 	}
// 	else {
// 		return x[]
// 	}
// }

// // Método Quase-Newton

// void QuasiNewlton(x1,x2,hess) {
// 	vector<float> grad = gradf(x1,x2);
// 	vector<float> x = x.push_back(x1).push_back(x2);
// 	vector<float> xk;
// 	matrix H = hess;
// 	if (grad.at(0) != 0 && grad.at(1) != 0) {
// 		vector<float> d = d.push_back(-(H.at(0)*grad.at(0)+H.at(1)*grad.at(1)).push_back(-(H.at(2)*grad.at(0)+H.at(3)*grad.at(1))
// 		//obter t>0
// 		xk.push_back(x1 + t*d.at(0)).push_back(x2 + t*d.at(1));
// 		// determinar proximo H
// 		vector<float> gradk = gradf(xk.at(0),xk.at(1));
// 		vector<float> p = p.push_back(xk.at(0) - x.at(0)).push_back(xk.at(1) - x.at(1));
// 		vector<float> q = q.push_back(gradk.at(0) - grad.at(0).push_back(gradk.at(1) - grad.at(1));
// 		matrix Hk;
// 		return QuasiNewlton(xk.at(0),xk.at(1));
// 	}	else {
// 		return x[];
// 	}
// }