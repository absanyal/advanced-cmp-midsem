/*
 * lib_harmonic.cpp
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
#include "common_globals.h"
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace boost::math;
using namespace std;

double alpha, omega;
double mass=1;
int number_of_mesh; double epsilon;
double offset = 1e-4;
double low_lim = -M_PI/2+offset;
double up_lim = M_PI/2-offset;
const int N4=10;

double fac(int n) {return factorial <double> (n);}
double psi(int i, double y);
double integrate_y(double (*func)(double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
double integrate_x(double y, double (*func_x)(double, double,int, int), int l, int n); //func_x= v_x
double v_x(double x, double y, int a, int b);
double v_d(double y, int l, int p, int n, int q);
double v_ex(double y, int l, int p, int n, int q);
string current_time_str(void);
void createfilename(string& filename, string option, int number_of_mesh, double omega);
void load_array_from_file(double direct[][N4][N4][N4],double exchange[][N4][N4][N4], ifstream& din, ifstream& ein);
void load_array_from_file(double matrixelems[][N4], ifstream& fin);
void generate_lhn_matrix(double omega, ofstream& fout);

double psi(int n, double x)
{
  double y = alpha*x;
  double result = 1/sqrt(pow(2,n)*fac(n))*sqrt(alpha)/pow(M_PI,0.25)*exp(-y*y/2)*hermite(n,y);
  return result;
}

double integrate_y(double (*func)(double, int, int, int, int), int l, int p, int n, int q ) //func = v_d/v_ex
{
  double trapez_sum;
  double fa, fb, v, step;
  int j;
  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func)(tan(low_lim),l,p,n,q)/(2.*pow(cos(low_lim),2)) ;
  fb=(*func)(tan(up_lim),l,p,n,q)/(2.*pow(cos(up_lim),2)) ;
  trapez_sum=0.;
  for (j=1; j <= number_of_mesh-1; j++)
  {
    v=j*step+low_lim;
    trapez_sum+=(*func)(tan(v),l,p,n,q)/pow(cos(v),2);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}


double v_d(double y, int l, int p, int n, int q)
{
  double integrated_x = integrate_x(y, &v_x, l, n);
  return integrated_x*psi(p,y)*psi(q,y);
}

double v_ex(double y, int l, int p, int n, int q)
{  cout << "Enter omega:";
  cin >> omega;
  alpha = sqrt(omega);
  double integrated_x = integrate_x(y, &v_x, l, q);
  return integrated_x*psi(p,y)*psi(n,y);
}

double integrate_x(double y, double (*func_x)(double, double,int, int), int l, int n) //func_x= v_x
{
  double trapez_sum;
  double fa, fb,u, step;
  int j;
  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func_x)(tan(low_lim),y,l,n)/(2.*pow(cos(low_lim),2));
  fb=(*func_x)(tan(up_lim),y,l,n)/(2.*pow(cos(up_lim),2));
  trapez_sum=0.;
  for (j=1; j <= number_of_mesh-1; j++)
  {
    u=j*step+low_lim;
    trapez_sum+=(*func_x)(tan(u),y,l,n)/pow(cos(u),2);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}


double v_x(double x, double y, int a, int b)
{
  return psi(a,x)*psi(b,x)/(epsilon+abs(x-y));
}


string current_time_str(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%Y-%m-%d-%I-%M-%S",timeinfo);
  string str(buffer);
  return str;
}

void createfilename(string& filename, string option, int number_of_mesh, double omega)
{
  filename ="data/"+option+"_harmonic(omega="+to_string(int(omega))+" ,mesh="+to_string(number_of_mesh)+","+current_time_str()+").txt";
}

void load_array_from_file(double direct[][N4][N4][N4],double exchange[][N4][N4][N4], ifstream& din, ifstream& ein)
{
  double d,ex;
  for(int l=0;l<N4;l++)
  {
    for(int p=0; p<N4; p++)
    {
      for(int n=0; n<N4; n++)
      {
        for(int q=0; q<N4; q++)
        {
          din >> d;
          ein >> ex;

          if(d<1e-4) direct[l][p][n][q]= 0.0;
          else  direct[l][p][n][q]=d;

          if(ex<1e-4)exchange[l][p][n][q]= 0.0;
          else  exchange[l][p][n][q]=ex;
        }
      }
    }
  }
}

void load_array_from_file(double matrixelems[][N4], ifstream& fin)
{
  double ld;
  for(int l=0; l<N4; l++)
  {
    for(int n=0; n<N4; n++)
    {
      fin >> ld;
      if(ld<1e-4)matrixelems[l][n]=0;
      else matrixelems[l][n]=ld;
    }
  }
}

double delta(int l, int n) {if(l==n) return 1.0; else return 0.0;}

void generate_lhn_matrix(double omega, ofstream& fout)
{
  for(int l=0; l<N4; l++)
  {
    for(int n=0; n<N4; n++)
      fout << (n+0.5)*omega*delta(l,n) << endl;
  }
}
