/*
 * lib_squarewell.cpp
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

using namespace std;

long double alpha, omega;
long double mass=1;
int number_of_mesh; long double epsilon;
long double low_lim = 0;
long double up_lim = 1;
const int N1=10;

long double psi(int i, long double y);
long double integrate_y(long double (*func)(long double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
long double integrate_x(long double y, long double (*func_x)(long double, long double,int, int), int l, int n); //func_x= v_x
long double v_x(long double x, long double y, int a, int b);
long double v_d(long double y, int l, int p, int n, int q);
long double v_ex(long double y, int l, int p, int n, int q);
string current_time_str(void);
void createfilename(string& filename, string option, int number_of_mesh, long double omega);
void load_array_from_file(long double direct[][N1][N1][N1],long double exchange[][N1][N1][N1], ifstream& din, ifstream& ein);
void load_array_from_file(long double matrixelems[][N1], ifstream& fin);
void generate_lhn_matrix(long double omega, ofstream& fout);

long double psi(int i, long double y)
{
  return sqrt(2)*sin((i+1)*M_PI*y);
}

long double integrate_y(long double (*func)(long double, int, int, int, int), int l, int p, int n, int q ) //func = v_d/v_ex
{
  long double trapez_sum;
  long double fa, fb, y, step;
  int j;
  step=(up_lim - low_lim)/((long double) number_of_mesh);
  fa=(*func)(low_lim,l,p,n,q);
  fb=(*func)(up_lim,l,p,n,q);
  trapez_sum=0.;
  for (j=1; j <= number_of_mesh-1; j++)
  {
    y=j*step+low_lim;
    trapez_sum+=(*func)(y,l,p,n,q);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

long double v_d(long double y, int l, int p, int n, int q)
{
  long double integrated_x = integrate_x(y, &v_x, l, n);
  return integrated_x*psi(p,y)*psi(q,y);
}

long double v_ex(long double y, int l, int p, int n, int q)
{
  long double integrated_x = integrate_x(y, &v_x, l, q);
  return integrated_x*psi(p,y)*psi(n,y);
}

long double integrate_x(long double y, long double (*func_x)(long double, long double,int, int), int l, int n) //func_x= v_x
{
  long double trapez_sum;
  long double fa, fb,x, step;
  int j;
  step=(up_lim - low_lim)/((long double) number_of_mesh);
  fa=(*func_x)((low_lim),y,l,n);
  fb=(*func_x)((up_lim),y,l,n);
  trapez_sum=0.;
  for (j=1; j <= number_of_mesh-1; j++)
  {
    x=j*step+low_lim;
    trapez_sum+=(*func_x)(x,y,l,n);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

long double v_x(long double x, long double y, int a, int b)
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

void createfilename(string& filename, string option, int number_of_mesh, long double omega)
{
  filename = "data/"+option+"_square_well(a=1,mesh="+to_string(number_of_mesh)+","+current_time_str()+").txt";
}

void load_array_from_file(long double direct[][N1][N1][N1],long double exchange[][N1][N1][N1], ifstream& din, ifstream& ein)
{
  long double d,ex;
  for(int l=0;l<N1;l++)
  {
    for(int p=0; p<N1; p++)
    {
      for(int n=0; n<N1; n++)
      {
        for(int q=0; q<N1; q++)
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

void load_array_from_file(long double matrixelems[][N1], ifstream& fin)
{
  long double ld;
  for(int l=0; l<N1; l++)
  {
    for(int n=0; n<N1; n++)
    {
      fin >> ld;
      if(ld<1e-4)matrixelems[l][n]=0;
      else matrixelems[l][n]=ld;
    }
  }
}

long double delta(int l, int n){if(l==n)return 1.0; else return 0.0;}

long double integrand(long double x, int l, int n) { return sin(l*M_PI*x)*sin(n*M_PI*x)*x*x; }

long double integrate_x2(long double (*func)(long double, int, int), int l, int n) //func = v_d/v_ex
{
  long double trapez_sum;
  long double fa, fb, y, step;
  int j;
  step=(up_lim - low_lim)/((long double) number_of_mesh);
  fa=(*func)(low_lim,l,n);
  fb=(*func)(up_lim,l,n);
  trapez_sum=0.0;
  for (j=1; j < number_of_mesh; j++)
  {
    y=j*step+low_lim;
    trapez_sum+=(*func)(y,l,n);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

void generate_lhn_matrix(long double omega, ofstream& fout)
{
  for(int i=0; i<N1; i++)
  {
    for(int j=0; j<N1; j++)
    {
      int l=i+1; int n=j+1;
      long double ld= 2*pow(omega,2)*integrate_x2(integrand,l,n)+ double(0.5*pow(n,2))*pow(M_PI,2)*delta(l,n);
      fout << ld << endl;
    }
  }
}
