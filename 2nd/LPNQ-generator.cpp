#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "common_globals.h"
#include "functions.h"

using namespace std;

ofstream dout,eout,fout;
const int N6=10;

int main()
{
  cout << "Enter the no of mesh pts: ";
  cin >> number_of_mesh;
  epsilon = 1/(double(number_of_mesh)*2.);
  string filename;

  cout << "Enter omega:";
  cin >> omega;
  alpha = sqrt(omega); //(In natural units with hbar=1.)

  int l,p,n,q;
  cin >> l >> p >> n >> q;
  cout << "Direct energy = " <<  integrate_y(&v_d,l, p, n,q) << endl;
  cout  << "Exchange energy = " <<  integrate_y(&v_ex,l, p, n,q) << endl;

  // createfilename(filename,"direct",number_of_mesh,omega); dout.open(filename);
  // createfilename(filename,"exchange",number_of_mesh,omega); eout.open(filename);
  // createfilename(filename,"matrixelems",number_of_mesh,omega); fout.open(filename);
  //
  // for(int l=0;l<N6;l++)
  // {
  //   for(int p=0; p<N6; p++)
  //   {
  //     for(int n=0; n<N6; n++)
  //     {
  //       for(int q=0; q<N6; q++)
  //       {
  //         dout <<  integrate_y(&v_d,l, p, n,q) << endl;
  //         eout <<  integrate_y(&v_ex,l, p, n,q) << endl;
  //       }
  //     }
  //    }
  //  }
  // generate_lhn_matrix(omega, fout);

   dout.close();
   eout.close();
   fout.close();

}
