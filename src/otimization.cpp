#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "otimization_functions.h"

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

float MathFunctions::phi(const vector<float> &X, const vector<float> &dX, const float t){
  vector<float> r = {X[0]+t*dX[0], X[1]+t*dX[1]};
  return funcObj(r);
}

vector<float> MathFunctions::gradF(const vector<float> &X) {
  vector<float> dX = vector<float>(2,0);
	dX[0] = (((pow(X[0],2))*(pow(X[1],2))-1)/((pow(X[0],3))*(pow(X[1],2))+X[0]));
	dX[1] = (((pow(X[0],2))*(pow(X[1],2))-1)/((pow(X[0],2))*(pow(X[1],3))+X[1]));
	return dX;
}