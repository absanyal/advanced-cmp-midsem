#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <lapacke.h>

using namespace std;
using namespace std::chrono;
using namespace Eigen;

typedef std::complex <double> cd;

#include "hermite_polynomial.hpp"

double Sqr(cd x){return (x*conj(x)).real();}
bool cEP(MatrixXcd A, VectorXcd& lambda, MatrixXcd& v)
{
  int N = A.cols();
  if (A.rows()!=N)  return false;
  v.resize(N,N);
  lambda.resize(N);

  int LDA = A.outerStride();
  int LDV = v.outerStride();
  int INFO = 0;
  cd* w = const_cast<cd*>(lambda.data());
  char Nchar = 'N';
  char Vchar = 'V';
  int LWORK = int(A.size())*4;
  VectorXcd WORK(LWORK);
  VectorXd RWORK(2*LDA);

  zgeev_(&Nchar, &Vchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, reinterpret_cast <__complex__ double*> (w), 0, &LDV, reinterpret_cast <__complex__ double*> (v.data()), &LDV,  reinterpret_cast <__complex__ double*> (WORK.data()), &LWORK, RWORK.data(), &INFO);

  for(int i=0; i<N; i++)
 	 v.col(i)=v.col(i)/v.col(i).unaryExpr(&Sqr).sum();

  return INFO==0;
}


int no_of_pts=1000;
const int number_of_mesh=100;
const double a = 1;
double low_lim = -8;
double up_lim = 8;
double dx = (up_lim - low_lim)/double(no_of_pts);
const double omega=1;
double alpha = 1/sqrt(omega);
const int no_of_sps = 10;
VectorXd point(no_of_pts+1);
MatrixXcd states(point.size(),no_of_sps);

void show_time(milliseconds begin_ms, milliseconds end_ms);
double fac(int n) {double prod=1.0; for(int i=n; i>1;i--) prod*=n; return prod;}
double V(double x) {return 0.5*pow(omega*x,2);}

double hermite(int n, double y)
{
    double x_vec[1];
    x_vec[0]=y;
    double* fx2_vec = h_polynomial_value ( 1, n, x_vec );
    double fx2 = fx2_vec[n];
    return fx2;
}


long double psi(int n, double x)
{
  long double y = alpha*x;
  long double result = 1/sqrt(pow(2,n)*fac(n))*sqrt(alpha)/pow(M_PI,0.25)*exp(-y*y/2)*hermite(n,y);
  return result;
}

bool compare(const pair<double, VectorXcd>&i, const pair<double, VectorXcd>&j) {return i.first < j.first;}

double filter(double x) {if(abs(x)<1e-3) return 0.0; else return x;}

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

double rho_H(double r)
{
  int n = int((r - low_lim)/dx);
  double rho=0.0;
  for(int i=0; i< no_of_sps; i++)
  {
    rho += (conj(states(n,i))*states(n,i)).real();
  }
  return rho;
}

double rho_HF(double r, double r_prime)
{
  double num=0.0; double denom=0.0;
  int n = (r - low_lim)/dx;
  int n_prime = (r_prime - low_lim)/dx;

  for(int i=0; i< no_of_sps; i++)
  {
    for(int j=0; j< no_of_sps; j++)
    {
      num += (conj(states(n,i))*states(n,i)*conj(states(n_prime,j))*states(n_prime,j)).real();
    }
    denom += (conj(states(n,i))*states(n,i)).real();
  }
  if(denom != 0) return num/denom; else return 0.0;
}

double integrand(double r, double r_prime)
{
   if(r!=r_prime) return (rho_H(r_prime)*rho_HF(r,r_prime))/abs(r - r_prime); else return 0.0;
 }

double integrate_rho(double r, double (*func_x)(double, double))
{
  double trapez_sum;
  double fa, fb,x, step;
  int j;
  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func_x)(r,low_lim);
  fb=(*func_x)(r,up_lim);
  trapez_sum=0.;
  for (j=1; j <= number_of_mesh-1; j++)
  {
    x=j*step+low_lim;
    trapez_sum+=(*func_x)(r,x);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

int main()
{
 	
  double maxdev;
  cout << "Enter tolerance: "; cin >> maxdev; 


  for(int i=0; i<= no_of_pts; i++) {point(i)=low_lim+i*dx;}
  for(int i=0; i< point.size(); i++)
  {
    for(int j=0; j<no_of_sps; j++)
      states(i,j) = psi(j,point(i));
  }

  for(int j=0; j<no_of_sps; j++)
      states.col(j) = states.col(j)/sqrt(states.col(j).unaryExpr(&Sqr).sum());

  cout << "Initial Normalization= " << states.col(1).unaryExpr(&Sqr).sum() << endl;

  MatrixXcd H = MatrixXcd::Zero(point.size(),point.size());
  for(int i=0; i<point.size(); i++)
  {
      int j = (i==point.size()-1)? 0 : i+1;
      H(i,j)= -1/(2*dx*dx);
      H(j,i)= -1/(2*dx*dx);
      H(i,i) = 1/(dx*dx)+ V(point(i));
  }


  VectorXcd v; MatrixXcd eigenvectors; VectorXd eigenvalues;
  int output_states = 2; int master_loop = 0;
  VectorXd oldeival= VectorXd::Zero(output_states);
  VectorXd neweival= VectorXd::Zero(output_states);
  ofstream fout("initialstate.txt");
   
   for(int i=0; i<point.size(); i++) fout << point(i) << " " << (states(i,0)).real() << " " << (states(i,1)).real() << endl;
   fout.close();


  for(; ; )
  {
    cout << "Loop-" << master_loop << "\n============================\n";
    milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    cEP(H,v,eigenvectors);
    milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    show_time(begin_ms, end_ms);

    eigenvalues = v.real();
    vector < pair<double,VectorXcd> > eigenspectrum;
    for(int i=0; i<point.size(); i++)
      eigenspectrum.push_back(make_pair(filter(eigenvalues(i)),filter(eigenvectors.col(i))));
    sort(eigenspectrum.begin(),eigenspectrum.end(),compare);
    eigenspectrum.resize(no_of_sps);

    for(int i=0; i<no_of_sps; i++) states.col(i)= eigenspectrum[i].second;   	
    for(int i=0; i<output_states; i++) neweival(i) = eigenspectrum[i].first;    	
    cout << "Eigenvalues are: " << neweival.transpose() << endl << endl;

    double max_deviation = (neweival - oldeival).cwiseAbs().maxCoeff();
    if(max_deviation < maxdev) break; else cout << "Max deviation = " << max_deviation << endl;

    for(int i=0; i<point.size(); i++) {H(i,i) = 1/(dx*dx) + V(point(i)) + integrate_rho(point(i),&integrand); cout << i << '\r';}
    for(int i=0; i<oldeival.size(); i++) oldeival(i)= eigenspectrum[i].first;
    master_loop++; cout << endl;
  }

  	  fout.open("finalstate.txt");
  	  for(int i=0; i<point.size(); i++) fout << point(i) << " " << (states(i,0)).real()  << " " << (states(i,1)).real() << endl;
  	  fout.close();

  	cout << "Normalization\n";

  	cout << sqrt(states.col(1).unaryExpr(&Sqr).sum()) << endl;

  	fout.open("state.txt");
  	fout << states.col(1) << endl;
  	fout.close();
}

void show_time(milliseconds begin_ms, milliseconds end_ms)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout << "Diagonalization took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << "Diagonalization took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << "Diagonalization took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << "Diagonalization took " << t << "time. Wrong t received.\n"; }
}
