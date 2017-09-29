#ifndef _FUNCTIONS_H_INCLUDED_
#define _FUNCTIONS_H_INCLUDED_

#include <iostream>
#include <cmath>
#include "common-globals.h"

using namespace std;

double psi(int i, double y);
double integrate_y(double (*func)(double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
double integrate_x(double y, double (*func_x)(double, double,int, int), int l, int n); //func_x= v_x
double v_x(double x, double y, int a, int b);
double v_d(double y, int l, int p, int n, int q);
double v_ex(double y, int l, int p, int n, int q);

#endif
