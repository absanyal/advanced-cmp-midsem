#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include "common-globals.h"
#include "functions.h"

using namespace std;
using namespace Eigen;

typedef std::complex <long double> cd;

// ifstream din("direct.txt");
// ifstream ein("exchange.txt");
ifstream din("direct_square_well.txt");
ifstream ein("exchange_square_well.txt");

const int N=3;
long double direct_energy[N][N][N][N];
long double exchange_energy[N][N][N][N];

cd LPNQ(int l, int n, MatrixXcd C)
{
  cd lpnq;
  for(int p=0; p<N; p++)
  {
    for(int q=0; q<N; q++)
    {
      cd c_sum=0;
      for(int mu=0; mu<N; mu++)  c_sum += conj(C(mu,p))*C(mu,q);
      lpnq += c_sum*(direct_energy[l][p][n][q]-exchange_energy[l][p][n][q]);
    }
  }
  return lpnq;
}

cd debug_csum(int l, int n, MatrixXcd C)
{
  cd lpnq;
  for(int p=0; p<N; p++)
  {
    for(int q=0; q<N; q++)
    {
      cd c_sum=0;
      for(int mu=0; mu<N; mu++)  c_sum += conj(C(mu,p))*C(mu,q);
      cout << c_sum << "\t" << direct_energy[l][p][n][q]-exchange_energy[l][p][n][q] << endl;
      lpnq+= c_sum;
    }
  }
  return lpnq;
}


double delta(int l, int n)
{
  if(l==n)
  return 1.0;
  else
   return 0.0;
}


int main()
{
  cout << "Enter omega: ";
  cin >> omega;
  //cout << "Run integration-generator for different omega input" << endl;

  load_array_from_file(direct_energy,exchange_energy,din,ein);

  MatrixXcd F = MatrixXd::Zero(N,N);
  MatrixXcd C = MatrixXcd::Identity(N,N);

for(int master_loop=1; master_loop<10; master_loop++)
{
  cout << "Loop-" << master_loop << ": " << endl << "----------------------\n";
  cout << C << endl << endl;

  for(int l=0; l<N; l++)
  {
    for(int n=0; n<N; n++)
    {
      F(l,n)=delta(l,n)*(double(n)+0.5)*omega+LPNQ(l,n,C);
    }
  }
  // cout << "debugging\n==============\n";
  // cd debug =  debug_csum(0,1,C);
  // cout << " debug_csum(0,1,C)= " << debug << "\n============" << '\n';

  ComplexEigenSolver <MatrixXcd> ces;
  ces.compute(F);
  C=ces.eigenvectors().transpose();
  // cout << ces.eigenvectors().inverse()*F*ces.eigenvectors() << endl;
  cout  << ces.eigenvalues().real().transpose() << endl;
}


}
