#ifndef _FUNCTIONS_H_INCLUDED_
#define _FUNCTIONS_H_INCLUDED_

#include <iostream>
#include <cmath>
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

#endif
