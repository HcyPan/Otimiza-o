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

  cout << endl << "Matrix A" << endl << endl;
  m_func.initializeMat(A);
  m_func.printMat(A);

  cout << endl << "Matrix B" << endl << endl;
  m_func.initializeMat(B);
  m_func.printMat(B);

  cout << endl << "Vector V" << endl << endl;
  m_func.initializeVec(V);
  m_func.printVec(V);

  cout << endl << "Matrix Multiplication Result" << endl << endl;
  matrix M = m_func.matMul(A,B);
  m_func.printMat(M);

  cout << endl <<  "Matrix Summed Result" << endl << endl;
  matrix S = m_func.matSum(A,B);
  m_func.printMat(S);

  cout << endl <<  "Matrix Summed Value Result" << endl << endl;
  matrix S2 = m_func.matSumValue(A,5);
  m_func.printMat(S2);

  cout << endl <<  "Matrix Vector Multiplication Result" << endl << endl;
  vector<float> VR = m_func.vecMatMul(V,A);
  m_func.printVec(VR);

  cout << endl <<  "Objective Function Result" << endl << endl;
  float objValue = m_func.funcObj(vector <float> {5,6});
  cout << objValue << endl;

  cout << endl <<  "Phi Function Result" << endl << endl;
  float phiValue = m_func.phi(vector <float> {5,6}, vector <float> {2,3}, 2);
  cout << phiValue << endl;

  cout << endl <<  "gradF Function Result" << endl << endl;
  vector<float> dX = m_func.gradF(vector <float> {5,6});
  m_func.printVec(dX);

  cout << endl <<  "Golden Section Function Result" << endl << endl;
  float goldSect = m_func.goldenSection(vector <float> {5,6}, vector <float> {2,3});
  cout << goldSect << endl;

  cout << endl <<  "Armijo Function Result" << endl << endl;
  float armijo = m_func.armijo(vector <float> {5,6}, vector <float> {2,3});
  cout << armijo << endl;

  cout << endl <<  "Gradient Function By Armijo Result" << endl << endl;
  vector<float> gradientArmijo = m_func.gradient(vector <float> {5,6}, MathFunctions::getStepFunction::Armijo);
  m_func.printVec(gradientArmijo);

  cout << endl <<  "Gradient Function By Golden Section Result" << endl << endl;
  vector<float> gradientGoldSec = m_func.gradient(vector <float> {5,6}, MathFunctions::getStepFunction::GoldenSect);
  m_func.printVec(gradientGoldSec);

  cout << endl <<  "Newton Function By Armijo Result" << endl << endl;
  vector<float> newtonArmijo = m_func.newton(vector <float> {5,6}, MathFunctions::getStepFunction::Armijo);
  m_func.printVec(newtonArmijo);

  cout << endl <<  "Newton Function By Golden Section Result" << endl << endl;
  vector<float> newtonGoldenSect = m_func.newton(vector <float> {5,6}, MathFunctions::getStepFunction::GoldenSect);
  m_func.printVec(newtonGoldenSect);

  matrix identity = m_func.initializeMatIdentity(2);

  cout << endl <<  "quasiNewton Function By Armijo Result" << endl << endl;
  vector<float> quasiNewtonArmijo = m_func.quasiNewton(vector <float> {5,6}, identity, MathFunctions::getStepFunction::Armijo);
  m_func.printVec(quasiNewtonArmijo);

  cout << endl <<  "quasiNewton Function By Golden Section Result" << endl << endl;
  vector<float> quasiNewtonGoldenSect = m_func.quasiNewton(vector <float> {5,6}, identity, MathFunctions::getStepFunction::GoldenSect);
  m_func.printVec(quasiNewtonGoldenSect);


  cout << endl << "Fnished" << endl;
}


int main(){
  allFunctions();
}