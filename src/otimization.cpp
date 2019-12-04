#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "otimization_functions.h"

#define GOLDENSECTIONERROR 1  
#define GOLDENSECTIONRO 1
#define ARMIJON 0.25
#define ARMIJOY 0.8


using namespace std;
using matrix = vector<vector<float> >;

matrix MathFunctions::matMul(const matrix &A, const matrix &B){

  const int n = A.size();     // a rows
  const int m = A[0].size();  // a cols
  const int p = B[0].size();  // b cols

  matrix C = matrix(n, vector<float>(p, 0));

  for (int i = 0; i < p; ++i){
      for (int j = 0; j < m; ++j){
          float aux = 0;
          for (int k = 0; k < n; ++k){
              aux += A[i][k] * B[k][j];
          }
          C[i][j] = aux;
      }
  }
  return C;
};

vector<float> MathFunctions::matVecMul(const vector<float> &A, const matrix &B){

  const int r = B.size();     // a rows
  const int c = B[0].size();  // a cols

  vector<float> V = vector<float>(c, 0);

  for (int i = 0; i < r; ++i){
      float aux = 0;
      for (int j = 0; j < c; ++j){
          aux += A[j] * B[i][j];
      }
      V[i] = aux;
  }
  return V;
};

float MathFunctions::vecMul(const vector<float> &A, const vector<float> &B){
  const int r = A.size();     // size of vector
  float aux = 0;

  for (int i = 0; i < r; ++i){
    aux += A[i]*B[i];
  }
  return aux;
}

void MathFunctions::vecMulFloat(vector<float> &A, float B){
  const int r = A.size();     // size of vector
  for (int i = 0; i < r; ++i){
    A[i] = A[i]*B;
  }
}

matrix MathFunctions::matSum(const matrix &A, const matrix &B){

  const int r = A.size();     // a rows
  const int c = A[0].size();  // a cols

  matrix S = matrix(r, vector<float>(c, 0));

  for (int i = 0; i < r; ++i){
      for (int j = 0; j < c; ++j){
          S[i][j] = A[i][j] + B[i][j];
      }
  }
  return S;
};

matrix MathFunctions::matSumValue(const matrix &A, const float B){

  const int r = A.size();     // a rows
  const int c = A[0].size();  // a cols

  matrix S = matrix(r, vector<float>(c, 0));

  for (int i = 0; i < r; ++i){
      for (int j = 0; j < c; ++j){
          S[i][j] = A[i][j] + B;
      }
  }
  return S;
};

void MathFunctions::initializeMat(matrix &A){
  const int c_start = rand() % 100 + 1;
  for (int i = 0; i < A.size(); i++){
    for (int j = 0; j < A[0].size(); j++){
        A[i][j] = i*j + c_start;
    }
  }
};

void MathFunctions::initializeVec(vector<float> &A){
  for (int i = 0; i < A.size(); i++){
      A[i] = rand() % 100 + 1;
  }
};

void MathFunctions::printMat(const matrix &A){
  for (int i = 0; i < A.size(); i++){
      for (int j = 0; j < A[0].size(); j++){
        cout << A[i][j] << " ";
      }
      cout << endl;
  }
};

void MathFunctions::printVec(const vector<float> &A){
  for (int i = 0; i < A.size(); i++){
        cout << A[i] << " ";
  }
  cout << endl;
};

float MathFunctions::funcObj(const vector<float> &X){
  float res = -log((X[1]*X[0])/(1+pow((X[0]*X[1]),2)));
 	return res;
}

float MathFunctions::phi(const vector<float> &X, const vector<float> &dir, const float t){
  vector<float> r = {X[0]+t*dir[0], X[1]+t*dir[1]};
  return funcObj(r);
}

vector<float> MathFunctions::gradF(const vector<float> &X) {
  vector<float> dX = vector<float>(2,0);
	dX[0] = (((pow(X[0],2))*(pow(X[1],2))-1)/((pow(X[0],3))*(pow(X[1],2))+X[0]));
	dX[1] = (((pow(X[0],2))*(pow(X[1],2))-1)/((pow(X[0],2))*(pow(X[1],3))+X[1]));
	return dX;
}

float MathFunctions::goldenSection(const vector<float> &X, const vector<float> &dir){
  float eps = GOLDENSECTIONERROR, ro = GOLDENSECTIONRO;
	float a = 0, s = ro , b = 2*ro;
  const float section = (3.0-sqrt(5.0))/2;
	vector<float> delta = {section, 1 - (section)} ;

  while (phi(X, dir, b) < phi(X, dir, s)){
    a = s;
    s = b;
    b = 2*b;
  }

  float ba = b-a;
  float u = a + delta[0]*ba;
  float v = a + delta[1]*ba;

  while ((b-a) > eps){
    if (phi(X, dir, u) < phi(X, dir, v)){
      b = v;
      v = u;
      u = a + delta[0]*(b-a);
    }else{
      a = u;
      u = v;
      v = a + delta[1] * (b-a);
    }
  }

  return (u+v) / 2;
}

float MathFunctions::armijo(const vector<float> &X, const vector<float> &dir){
	float Y = ARMIJOY, N = ARMIJON;

  float t = 1;

	vector<float> grad = gradF(X);

	while (phi(X, dir, t) > funcObj(X) + N*t*(vecMul(grad, dir))){
		t = Y*t;
	}

  return t;
}

vector<float> MathFunctions::gradient(const vector<float> &X, int stepFunction){
  vector<float> grad = gradF(X);

  float t = 0;

  if (grad[0] != 0 && grad[1] != 0){
    vecMulFloat(grad, -1);
    if (stepFunction == 1){
      t = armijo(X, grad);
    }else{
      t = goldenSection(X, grad);
    }
    vecMulFloat(grad, t);

    return gradient(vector<float> {X[0] + grad[0], X[1] + grad[1]}, stepFunction);
  }
  return X;
}

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