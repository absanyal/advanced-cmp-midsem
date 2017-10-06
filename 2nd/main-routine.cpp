#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <sstream>
#include "common_globals.h"
#include "functions.h"

using namespace std;
using namespace Eigen;

typedef std::complex <long double> cd;

const int Nb=10;
const int N=2;

long double direct_energy[Nb][Nb][Nb][Nb];
long double exchange_energy[Nb][Nb][Nb][Nb];
long double matrixelems[Nb][Nb];

cd LPNQ(int l, int n, MatrixXcd C)
{
  cd lpnq;
  for(int p=0; p<Nb; p++)
  {
    for(int q=0; q<Nb; q++)
    {
      cd c_sum=0;
      for(int mu=0; mu<N; mu++)  c_sum += conj(C(mu,p))*C(mu,q);
      lpnq += c_sum*(2*direct_energy[l][p][n][q]-exchange_energy[l][p][n][q]);
    }
  }
  return lpnq;
}

bool compare(const pair<double, VectorXcd>&i, const pair<double, VectorXcd>&j) {return i.first < j.first;}

double filter(double x) {if(x<1e-3) return 0.0; else return x;}

VectorXcd filter(VectorXcd v)
{
  VectorXcd v_filtered(v.size());
  for(int i=0; i<v.size(); i++)
  {
    v_filtered(i).real(filter(v(i).real()));
    v_filtered(i).imag(filter(v(i).imag()));
  }
  return v_filtered;
}

int main(int argc, char* argv[])
{

  if(argc !=2) {cout << "pass proper arguments to main()\n"; exit(1);}
  istringstream ss(argv[1]); int choice;
  if (!(ss >> choice)) cerr << "Invalid number " << argv[1] << '\n';

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

  MatrixXcd F = MatrixXcd::Zero(Nb,Nb);
  MatrixXcd C = MatrixXcd::Zero(N,Nb);
  MatrixXcd Cn= MatrixXcd::Zero(N,Nb);
  int number_of_loops;
  cout << "Enter the number_of_loops: "; cin >> number_of_loops;

for(int master_loop=1; master_loop<number_of_loops; master_loop++)
{
  cout << "Loop-" << master_loop << ": " << endl << "----------------------\n";

  for(int l=0; l<Nb; l++)
  {
    for(int n=0; n<Nb; n++)
    {
      F(l,n)=matrixelems[l][n] + LPNQ(l,n,C);
    }
  }

  ComplexEigenSolver <MatrixXcd> ces; ces.compute(F);
  vector < pair<double,VectorXcd> > eigenspectrum;

  for(int i=0; i<Nb; i++)
    eigenspectrum.push_back(make_pair(filter(ces.eigenvalues().real()[i]),filter(ces.eigenvectors().col(i))));

  sort(eigenspectrum.begin(),eigenspectrum.end(),compare);
  eigenspectrum.resize(N);

  for(int i=0; i<N; i++)
  {
    C.row(i) = eigenspectrum[i].second.transpose()*0.75+ C.row(i)*0.25;
    cout << eigenspectrum[i].first << " \t ";
  }
  cout << endl;
  cout << C.real() << endl << endl;
}

  cout << "Final Fock Matrix: \n";

  cout << F.real() << endl << endl;
  din.close();
  ein.close();
  fin.close();

}
