#include <iostream>
#include <cmath>
#include "common.h"

using namespace std;

double psi(int i, double y);
double integrate_y(double (*func)(double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
double integrate_x(double y, double (*func_x)(double, double,int, int), int l, int n); //func_x= v_x
double v_x(double x, double y, int a, int b);
double v_d(double y, int l, int p, int n, int q);
double v_ex(double y, int l, int p, int n, int q);


double psi(int i, double y)
{
  double x = sqrt(alpha)*y;
  double result= pow((alpha/M_PI),0.25)*exp(-x*x/2);
  if(i==0) return result*1;
  else if(i==1) return result*sqrt(2)*x;
  else if(i==2) return result*(2*x*x-1)/sqrt(2);
  else { cout << "wrong psi called" << endl; exit(1); }
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
{
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
