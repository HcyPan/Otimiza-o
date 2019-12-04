#include "otimization_functions.h"
#include <vector>
#include <iostream> 

using namespace std;
using matrix = vector<vector<float> >;

void allFunctions(){
  cout << "Starting" << endl;
  MathFunctions m_func;
  matrix A = matrix(5, vector<float>(5, 0));
  matrix B = matrix(5, vector<float>(5, 0));
  vector<float> V = vector<float>(5, 0);

  m_func.initializeMat(A);
  cout << endl << "Matrix A" << endl << endl;
  m_func.printMat(A);

  m_func.initializeMat(B);
  cout << endl << "Matrix B" << endl << endl;
  m_func.printMat(B);

  m_func.initializeVec(V);
  cout << endl << "Vector V" << endl << endl;
  m_func.printVec(V);

  matrix M = m_func.matMul(A,B);
  cout << endl << "Matrix Multiplication Result" << endl << endl;
  m_func.printMat(M);

  matrix S = m_func.matSum(A,B);
  cout << endl <<  "Matrix Summed Result" << endl << endl;
  m_func.printMat(S);

  matrix S2 = m_func.matSumValue(A,5);
  cout << endl <<  "Matrix Summed Value Result" << endl << endl;
  m_func.printMat(S2);

  vector<float> VR = m_func.matVecMul(V,A);
  cout << endl <<  "Matrix Vector Multiplication Result" << endl << endl;
  m_func.printVec(VR);

  float objValue = m_func.funcObj(vector <float> {5,6});
  cout << endl <<  "Objective Function Result" << endl << endl;
  cout << objValue << endl;

  float phiValue = m_func.phi(vector <float> {5,6}, vector <float> {2,3}, 2);
  cout << endl <<  "Phi Function Result" << endl << endl;
  cout << phiValue << endl;

  vector<float> dX = m_func.gradF(vector <float> {5,6});
  cout << endl <<  "gradF Function Result" << endl << endl;
  m_func.printVec(dX);

  cout << endl << "Fnished" << endl;
}


int main(){
  allFunctions();
}