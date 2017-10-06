#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const int N = 100;
const double a = 1;
double low_lim = -2;
double up_lim = 2;
double dx = (up_lim - low_lim)/double(N);

double V(double x) 

int main()
{
  double point[N];
  for(int i=0; i<N; i++) {point[i]=low_lim+i*(up_lim - low_lim)/double(N)}

  MatrixXcd H = MatrixXcd::Zero(N,N);
  auto it = v.begin();
  for(int i=0; i<N; i++)
  {
      H(i,i+1)= -1/(2*dx*dx);
      H(i+1,i)= -1/(2*dx*dx);
      H(i,i) = 1/(dx*dx) + V(point[i]);
  }


}
