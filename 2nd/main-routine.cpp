#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include "common_globals.h"
#include "functions.h"

using namespace std;
using namespace Eigen;

typedef std::complex <long double> cd;

const int N3=10;

long double direct_energy[N3][N3][N3][N3];
long double exchange_energy[N3][N3][N3][N3];
long double matrixelems[N3][N3];

cd LPNQ(int l, int n, MatrixXcd C)
{
  cd lpnq;
  for(int p=0; p<N3; p++)
  {
    for(int q=0; q<N3; q++)
    {
      cd c_sum=0;
      for(int mu=0; mu<N3; mu++)  c_sum += conj(C(mu,p))*C(mu,q);
      lpnq += c_sum*(2*direct_energy[l][p][n][q]-exchange_energy[l][p][n][q]);
    }
  }
  return lpnq;
}


int main()
{
  int choice;
  cout << "Start with Harmonic Basis States(1) or Square Well Basis States(2): ";
  cin >> choice;

  ifstream din, ein, fin;

  if(choice==1)
  {
    din.open("direct_harmonic.txt");
    ein.open("exchange_harmonic.txt");
    fin.open("matrixelems_harmonic.txt");
  }
  else if(choice==2)
  {
    din.open("direct_square_well.txt");
    ein.open("exchange_square_well.txt");
    fin.open("matrixelems_square_well.txt");
  }
  else
  {
    cout << "Come on! Be professional and stop trolling the program" << endl;
    exit(1);
  }

  load_array_from_file(direct_energy,exchange_energy,din,ein);
  load_array_from_file(matrixelems,fin);

  MatrixXcd F = MatrixXcd::Zero(N3,N3);
  MatrixXcd C = MatrixXcd::Identity(N3,N3);

  for(int master_loop=1; master_loop<10; master_loop++)
  {
    cout << "Loop-" << master_loop << ": " << endl << "----------------------\n";
    cout << C << endl << endl;

    for(int l=0; l<N3; l++)
    {
      for(int n=0; n<N3; n++)
      {
        F(l,n)=matrixelems[l][n] + LPNQ(l,n,C);
      }
    }
    // if(master_loop==1)
    //  cout << "Initial F matrix is:\n"<< F << endl << endl << "and F(0,0) addends are: " << matrixelems[0][0] << " " << LPNQ(0,0,C) << endl;

    ComplexEigenSolver <MatrixXcd> ces;
    ces.compute(F);
    C=ces.eigenvectors().transpose();
    cout  << ces.eigenvalues().real().transpose() << endl;

  }

  din.close();
  ein.close();
  fin.close();

}
