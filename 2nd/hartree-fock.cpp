#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace boost::math;

const int N = 500;
const int number_of_mesh=100;
const double a = 1;
double low_lim = -4;
double up_lim = 4;
double dx = (up_lim - low_lim)/double(N);
const double omega=1;
double alpha = 1/sqrt(omega);
const int no_of_states = 10;
double h = (up_lim - low_lim)/double(N-1);

VectorXd point(N); MatrixXcd states(N,no_of_states);

void show_time(milliseconds begin_ms, milliseconds end_ms);
double fac(int n) {return factorial <double> (n);}
double V(double x) {return pow(omega*x,2);}

long double psi(int n, double x)
{
  long double y = alpha*x;
  long double result = 1/sqrt(pow(2,n)*fac(n))*sqrt(alpha)/pow(M_PI,0.25)*exp(-y*y/2)*hermite(n,y);
  return result;
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

double rho_H(double r)
{
  int n = (r - low_lim)/h;
  double rho=0.0;
  for(int i=0; i< no_of_states; i++)
  {
    rho += (conj(states(n,i))*states(n,i)).real();
  }
  return rho;
}

double rho_HF(double r, double r_prime)
{
  double num=0.0; double denom=0.0;
  int n = (r - low_lim)/h;
  int n_prime = (r_prime - low_lim)/h;

  for(int i=0; i< no_of_states; i++)
  {
    for(int j=0; j< no_of_states; j++)
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
  int no_of_loops;
  cout << "Enter the no_of_loops: ";
  cin >> no_of_loops;

  for(int i=0; i<N; i++) {point(i)=low_lim+i*h;}

  for(int i=0; i<N; i++)
  {
    for(int j=0; j<no_of_states; j++)
      states(i,j) = psi(j,point(i));
  }


  MatrixXcd H = MatrixXcd::Zero(N,N);

  for(int i=0; i<N; i++)
  {
      cout.flush();
      int j = (i==N-1)? 0 : i+1;
      cout << i << " " << j << endl;
      H(i,j)= -1/(2*dx*dx);
      H(j,i)= -1/(2*dx*dx);
      H(i,i) = 1/(dx*dx)+ V(point(i)) + integrate_rho(point(i),&integrand);
  }

  ComplexEigenSolver <MatrixXcd> ces;

  for(int master_loop=0; master_loop<no_of_loops; master_loop++)
  {
    cout << "Loop-" << master_loop << "\n============================\n";
    milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    ces.compute(H);
    milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    show_time(begin_ms, end_ms);

    vector < pair<double,VectorXcd> > eigenspectrum;
    for(int i=0; i<N; i++)
      eigenspectrum.push_back(make_pair(filter(ces.eigenvalues().real()[i]),filter(ces.eigenvectors().col(i))));

    sort(eigenspectrum.begin(),eigenspectrum.end(),compare);
    eigenspectrum.resize(no_of_states);
    for(int i=0; i<no_of_states; i++) states.col(i)= eigenspectrum[i].second;

    for(int i=0; i<no_of_states; i++) H(i,i) = 1/(dx*dx) + V(point(i)) + integrate_rho(point(i),&integrand);
    cout << ces.eigenvalues()[0].real() << "\t" << ces.eigenvalues()[1].real() << endl;
  }
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
