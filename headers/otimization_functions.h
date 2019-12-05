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

    bool isVecEqual(const vector<float> &A, vector<float> &B);

    matrix matMul(const matrix &A, const matrix &B);
    vector<float> vecMatMul(const vector<float> &A, const matrix &B);
    vector<float> matVecMul(const matrix &A, const vector<float> &B);
    float vecMul(const vector<float> &A, const vector<float> &B);
    vector<float> vecMulFloat(vector<float> &A, float B);
    vector<float> vecMagicOperation(const vector<float> &A, const vector<float> &B, float C);

    matrix hessianInverted(const vector<float> &X);

    matrix matSum(const matrix &A, const matrix &B);
    matrix matSumValue(const matrix &A, const float B);
    
    float funcObj(const vector<float> &X);
    float phi(const vector<float> &X, const vector<float> &dir, float t);
    vector<float> gradF(const vector<float> &X);
    float goldenSection(const vector<float> &X,  const vector<float> &dir);
    float armijo(const vector<float> &X,  const vector<float> &dir);
    vector<float> gradient(const vector<float> &X, int stepFunction, int iteration = 0);
    vector<float> newton(const vector<float> &X, int stepFunction, int iteration = 0);
};

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