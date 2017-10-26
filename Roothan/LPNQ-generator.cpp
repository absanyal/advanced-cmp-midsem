/*
 * LPNQ-generator.cpp
 *
 * Copyright 2017 Bineet Dash <bineet@bineet-ubuntu>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

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
