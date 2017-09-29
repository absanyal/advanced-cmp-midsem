#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "common-globals.h"
#include "functions.h"

using namespace std;
ofstream dout;
ofstream eout;

int main()
{
  cout << "Enter the no of mesh pts: ";
  cin >> number_of_mesh;
  epsilon = 1/(float(number_of_mesh)*2.);

  cout << "Enter omega:";
  cin >> omega;
  alpha = sqrt(omega); //In natural units with hbar=1.

  // cout <<  integrate_y(&v_d,0,1,2,1) << endl;
  // cout <<  integrate_y(&v_ex,0,1,2,1) << endl;

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
     }
   }

}
