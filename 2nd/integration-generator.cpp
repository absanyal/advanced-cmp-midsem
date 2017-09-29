#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "common-globals.h"
#include "functions.h"

using namespace std;
ofstream dout;
ofstream eout;

long double alpha, omega;
long double mass=1;
int number_of_mesh; long double epsilon;
long double offset = 1e-4;
long double low_lim = -M_PI/2+offset;
long double up_lim = M_PI/2-offset;

int main()
{
  cout << "Enter the no of mesh pts: ";
  cin >> number_of_mesh;
  epsilon = 1/(float(number_of_mesh)*2.);

  cout << "Enter omega:";
  cin >> omega;
  alpha = sqrt(omega); //In natural units with hbar=1.

  string filename;
  filename = "direct_omega="+to_string(omega)+"mesh="+to_string(number_of_mesh)+current_time_str(); dout.open(filename);
  filename = "exchange_omega="+to_string(omega)+"mesh="+to_string(number_of_mesh)+current_time_str(); eout.open(filename);

  for(int l=0;l<3;l++)
  {
    for(int p=0; p<3; p++)
    {
      for(int n=0; n<3; n++)
      {
        for(int q=0; q<3; q++)
        {
          dout <<  integrate_y(&v_d,l, p, n,q) << endl;
          eout <<  integrate_y(&v_ex,l, p, n,q) << endl;
        }
      }
      cout << "p=" << p << "is done!";
    }
    cout << "l=" << l << "is done!";
  }

}
