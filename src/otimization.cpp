#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "otimization_functions.h"

#define GOLDENSECTIONERROR 1  
#define GOLDENSECTIONRO 1
#define ARMIJON 0.25
#define ARMIJOY 0.8
#define RANDOMMAX 1000
#define MAXITERS 20000

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
}

vector<float> MathFunctions::vecMatMul(const vector<float> &A, const matrix &B){

  const int r = B.size();     // a rows
  const int c = B[0].size();  // a cols

  vector<float> V = vector<float>(c, 0);

  for (int i = 0; i < r; ++i){
      float aux = 0;
      for (int j = 0; j < c; ++j){
          aux += A[j] * B[j][i];
      }
      V[i] = aux;
  }
  return V;
}

vector<float> MathFunctions::matVecMul(const matrix &A, const vector<float> &B){

  const int r = A.size();     // a rows
  const int c = A[0].size();  // a cols

  vector<float> V = vector<float>(c, 0);

  for (int i = 0; i < r; ++i){
      float aux = 0;
      for (int j = 0; j < c; ++j){
          aux += B[j] * A[i][j];
      }
      V[i] = aux;
  }
  return V;
}

float MathFunctions::vecMul(const vector<float> &A, const vector<float> &B){
  const int r = A.size();     // size of vector
  float aux = 0;

  for (int i = 0; i < r; ++i){
    aux += A[i]*B[i];
  }
  return aux;
}

vector<float> MathFunctions::vecMulFloat(vector<float> &A, float B){
  const int r = A.size();     // size of vector
  vector<float> V = vector<float>(r, 0);

  for (int i = 0; i < r; ++i){
    V[i] = A[i]*B;
  }

  return V;
}

vector<float> MathFunctions::vecMagicOperation(const vector<float> &A, const vector<float> &B, float C){
  const int r = A.size();     // size of vector
  
  vector<float> V = vector<float>(r, 0);

  for (int i = 0; i < r; ++i){
    V[i] = A[i] + C*B[i];
  }

  return V;
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
}

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
}

void MathFunctions::initializeMat(matrix &A){
  for (int i = 0; i < A.size(); i++){
    for (int j = 0; j < A[0].size(); j++){
        A[i][j] = rand() % RANDOMMAX + 1;
    }
  }
}

void MathFunctions::initializeVec(vector<float> &A){
  for (int i = 0; i < A.size(); i++){
      A[i] = rand() % RANDOMMAX + 1;
  }
}

void MathFunctions::printMat(const matrix &A){
  for (int i = 0; i < A.size(); i++){
      for (int j = 0; j < A[0].size(); j++){
        cout << A[i][j] << " ";
      }
      cout << endl;
  }
}

void MathFunctions::printVec(const vector<float> &A){
  for (int i = 0; i < A.size(); i++){
        cout << A[i] << " ";
  }
  cout << endl;
}

float MathFunctions::funcObj(const vector<float> &X){
  float res = -log((X[1]*X[0])/(1+pow((X[0]*X[1]),2)));
 	return res;
}

float MathFunctions::phi(const vector<float> &X, const vector<float> &dir, const float t){
  vector<float> r = vecMagicOperation(X, dir, t);
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
	vector<float> delta = {section, 1 - (section)};

  while (phi(X, dir, b) < phi(X, dir, s)){
    a = s;
    s = b;
    b = 2*b;
  }

  float ba = b-a;

  vector<float> newDelta = vecMulFloat(delta, ba);

  float u = a + newDelta[0];
  float v = a + newDelta[1];

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

bool MathFunctions::isVecEqual(const vector<float> &A, vector<float> &B){
  const int r = A.size();  //size of vec
  bool C = true;

  for (int i = 0; i < r; i++){
      C = C && (A[i] == B[i]);
  }

  return C;
}

vector<float> MathFunctions::gradient(const vector<float> &X, int stepFunction, int iteration){
  iteration++;

  if (iteration == MAXITERS){
    cout << "Finished after max iteration steps: " << iteration << endl;
    return X;
  }

  vector<float> grad = gradF(X);

  float t = 0;

  if (grad[0] != 0 && grad[1] != 0){
    vector<float> dir = vecMulFloat(grad, -1);
    if (stepFunction == 1){
      t = armijo(X, dir);
    }else{
      t = goldenSection(X, dir);
    }
    
    vector<float> returnValue = vecMagicOperation(X, dir, t);

    if (isVecEqual(X, returnValue)){ //checks if new X is equal
      cout << "Finished after X equals to new X. " << endl;
      return X;
    }

    return gradient(returnValue, stepFunction, iteration);
  }
  cout << "Finished because grad is 0" << endl;
  return X;
}

matrix MathFunctions::hessianInverted(const vector<float> &X){
  float denominator = pow(X[0],6) * pow(X[1],6)  - 9 * pow(X[0],4) * pow(X[1],4) + 7 * pow(X[0],2) * pow(X[1],2) + 1;
  float secondDiogonal = (4 * pow(X[0],3) * pow(X[1],3) * (pow(X[0],2) * pow(X[1],2) + 1))/denominator;
  float firstDiagonal11 = (-pow(X[1],6) * pow(X[0],8) + 3 * pow(X[1],4) * pow(X[0], 6) + 5* pow(X[1], 2) * pow(X[0], 4) + pow(X[0],2))/denominator;
  float firstDiagonal22 = (-pow(X[0],6) * pow(X[1],8) + 3 * pow(X[0],4) * pow(X[1], 6) + 5* pow(X[0], 2) * pow(X[1], 4) + pow(X[1],2))/denominator;

  return {{firstDiagonal11, secondDiogonal}, {secondDiogonal, firstDiagonal22}};
}


vector<float> MathFunctions::newton(const vector<float> &X, int stepFunction, int iteration){
  iteration++;
  
  if (iteration == MAXITERS){
    cout << "Finished after max iteration steps: " << iteration << endl;
    return X;
  }

  vector<float> grad = gradF(X);

  float t = 0;

  if (grad[0] != 0 && grad[1] != 0){
    matrix HI = hessianInverted(X);
    vector<float> HIgrad = matVecMul(HI, grad);
    vector<float> dir = vecMulFloat(HIgrad, -1);
    if (stepFunction == 1){
      t = armijo(X, dir);
    }else{
      t = goldenSection(X, dir);
    }
    
    vector<float> returnValue = vecMagicOperation(X, dir, t);

    if (isVecEqual(X, returnValue)){ //checks if new X is equal
      cout << "Finished after X equals to new X" << endl;
      return X;
    }

    return newton(returnValue, stepFunction, iteration);
  }
  cout << "Finished because grad is 0" << endl;
  return X;
}

