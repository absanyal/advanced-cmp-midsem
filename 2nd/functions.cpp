#include <iostream>
#include <cmath>
#include <fstream>
#include "common-globals.h"

using namespace std;

long double psi(int i, long double y);
long double integrate_y(long double (*func)(long double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
long double integrate_x(long double y, long double (*func_x)(long double, long double,int, int), int l, int n); //func_x= v_x
long double v_x(long double x, long double y, int a, int b);
long double v_d(long double y, int l, int p, int n, int q);
long double v_ex(long double y, int l, int p, int n, int q);
string current_time_str(void);
void load_array_from_file(long double**** direct,long double**** exchange, ifstream& din, ifstream& ein);

long double psi(int i, long double y)
{
  long double x = sqrt(alpha)*y;
  long double result= pow((alpha/M_PI),0.25)*exp(-x*x/2);
  if(i==0) return result*1;
  else if(i==1) return result*sqrt(2)*x;
  else if(i==2) return result*(2*x*x-1)/sqrt(2);
  else { cout << "wrong psi called" << endl; exit(1); }
}

long double integrate_y(long double (*func)(long double, int, int, int, int), int l, int p, int n, int q ) //func = v_d/v_ex
{
  long double trapez_sum;
  long double fa, fb, v, step;
  int j;
  step=(up_lim - low_lim)/((long double) number_of_mesh);
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
  long double fa, fb,u, step;
  int j;
  step=(up_lim - low_lim)/((long double) number_of_mesh);
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
  strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
  string str(buffer);
  return str;
}

void load_array_from_file(long double**** direct,long double**** exchange, ifstream& din, ifstream& ein)
{
  for(int l=0;l<3;l++)
  {
    for(int p=0; p<3; p++)
    {
      for(int n=0; n<3; n++)
      {
        for(int q=0; q<3; q++)
        {
          din >> direct[l][p][n][q];
          ein >> exchange[l][p][n][q];
        }
      }
      cout << "p=" << p << "is done!";
    }
    cout << "l=" << l << "is done!";
  }
}
