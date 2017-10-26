/*
 * direct-exchange.cpp
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
#include <fstream>
#include <cmath>
#include <chrono>
#include "functions.h"
#include "common_globals.h"

using namespace std;
using namespace std::chrono;

double phi(double x, double y){return 1/sqrt(2)*(psi(0,x)*psi(1,y)-psi(0,y)*psi(1,x)); }

double first_order_correction(double x, double y) { return pow(phi(x,y),2)/(abs(x-y)+epsilon);}

double integrate_x_perturbation(double y, double (*func_x)(double, double)) //func_x=first_order_correction
{
  double trapez_sum;
  double fa, fb,u, step;

  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func_x)(tan(low_lim),y)/(2.*pow(cos(low_lim),2));
  fb=(*func_x)(tan(up_lim),y)/(2.*pow(cos(up_lim),2));
  trapez_sum=0.;
  for (int j=1; j < number_of_mesh; j++)
  {
    u=j*step+low_lim;
    trapez_sum+=(*func_x)(tan(u),y)/pow(cos(u),2);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

double integrand_y_perturbation(double y)
{
  return integrate_x_perturbation(y,&first_order_correction);
}

double integrate_y_perturbation(double (*func)(double)) //func = v_d/v_ex
{
  double trapez_sum;
  double fa, fb, v, step;

  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func)(tan(low_lim))/(2.*pow(cos(low_lim),2)) ;
  fb=(*func)(tan(up_lim))/(2.*pow(cos(up_lim),2)) ;
  trapez_sum=0.;

  for (int j=1; j < number_of_mesh; j++)
  {
    v=j*step+low_lim;
    trapez_sum+=(*func)(tan(v))/pow(cos(v),2);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

int main()
{
  number_of_mesh=500;
  epsilon = 1/(double(number_of_mesh)*2.);
  //
  // cout << "Enter omega:";
  // cin >> omega;

  // ofstream fout("perturbtion_data.txt");

  // for(omega=0.01; omega <3.0; omega+= 0.01)
  // {
  //   alpha = sqrt(omega);
  //   double res= integrate_y_perturbation(&integrand_y_perturbation);
  //   fout << omega << " " << res << endl;
  //   cout << omega << endl;
  // }



}
