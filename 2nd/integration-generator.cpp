#include <iostream>
#include <cmath>
#include "common-globals.h"
#include "functions.h"

using namespace std;

double alpha, omega, mass;
int number_of_mesh; double epsilon;
double offset = 1e-4;
double low_lim = -M_PI/2+offset;
double up_lim = M_PI/2-offset;

int main()
{
  cout << "Enter the no of mesh pts: ";
  cin >> number_of_mesh;
  epsilon = 1/(float(number_of_mesh)*2.);

  cout << "Enter omega and mass: ";
  cin >> omega >> mass;
  alpha = sqrt(mass*omega); //In natural units with hbar=1.

  int l,p,n,q;
  cout << "Enter l,p,n,q: ";
  cin >> l >> p >> n >>  q;
  // cout << "epsilon= " << epsilon << " \n" << "Direct Energy= " << integrate_y(&v_d,l, p, n,q) << endl;
  // cout << "Exchange Energy= " << integrate_y(&v_ex, l,p ,n,q) << endl;
  cout << "Diff =" <<  integrate_y(&v_d,l, p, n,q)- integrate_y(&v_ex,l, p, n,q) << endl;
}
