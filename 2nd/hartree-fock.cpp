#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace Eigen;
using namespace boost::math;

const int N = 100;
const number_of_mesh=100;
const double a = 1;
double low_lim = -2;
double up_lim = 2;
double dx = (up_lim - low_lim)/double(N);
const double omega=1;
double alpha = 1/sqrt(omega);
const int no_of_states = 10;
VectorXd point(N); MatrixXcd states(N,N);


double fac(int n) {return factorial <double> (n);}
double V(double x) {return pow(omega*x,2);}

long double psi(int n, double x)
{
  long double y = alpha*x;
  long double result = 1/sqrt(pow(2,n)*fac(n))*sqrt(alpha)/pow(M_PI,0.25)*exp(-y*y/2)*hermite(n,y);
  return result;
}

double rho_H(double r)
{
  double rho=0.0;
  for(int i=0; i<= no_of_states; i++)
  {
    rho += conj(states(i,r))*states(i,r).real();
  }
  return rho;
}

double rho_HF(double r, double r_prime)
{
  double num=0.0; double denom=0.0;
  for(int i=0; i<= no_of_states; i++)
  {
    for(int j=0; j<= number_of_loops; j++)
    {
      num += conj(states(i,r))*states(i,r)*conj(states(j,r_prime))*states(j,r_prime).real();
    }
    denom += conj(states(i,r))*states(i,r).real();
  }
  (denom != 0): return num/denom ? return 0.0;
}

double integrand(double r, double r_prime)
{
   (r==r_prime): return (rho_H(r_prime)*rho_HF(r,r_prime))/abs(r - r_prime) ? return 0;
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

  for(int i=0; i<N; i++) {point(i)=low_lim+i*(up_lim - low_lim)/double(N);}
  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
      states(i,j) = psi(j,point(i));
  }

  MatrixXcd H = MatrixXcd::Zero(N,N);

  for(int i=0; i<N; i++)
  {
      H(i,i+1)= -1/(2*dx*dx);
      H(i+1,i)= -1/(2*dx*dx);
      H(i,i) = 1/(dx*dx) + V(point(i)) + integrate_rho(point(i),&integrand);
  }

  ComplexEigenSolver <MatrixXcd> ces;

  for(int master_loop=1; master_loop<no_of_loops; master_loop++)
  {
    ces.compute(H);
    for(int i=0; i<no_of_states; i++) states.col(i)=ces.eigenvectors.col(i);
    for(int i=0; i<no_of_states; i++) H(i,i) = 1/(dx*dx) + V(point(i)) + integrate_rho(point(i),&integrand);
    cout << ces.eigenvalues(0) << "\t" << ces.eigenvalues(0) << endl;
  }




}
