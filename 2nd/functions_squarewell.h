#ifndef _FUNCTIONS_SQUAREWELL_H_INCLUDED_
#define _FUNCTIONS_SQUAREWELL_H_INCLUDED_

#include <iostream>
#include <cmath>
#include "common_globals.h"

using namespace std;
const int N2=3;

long double psi(int i, long double y);
long double integrate_y(long double (*func)(long double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
long double integrate_x(long double y, long double (*func_x)(long double, long double,int, int), int l, int n); //func_x= v_x
long double v_x(long double x, long double y, int a, int b);
long double v_d(long double y, int l, int p, int n, int q);
long double v_ex(long double y, int l, int p, int n, int q);
string current_time_str(void);
void createfilename(string& filename, string option, int number_of_mesh, long double omega);
void load_array_from_file(long double direct[][N2][N2][N2],long double exchange[][N2][N2][N2], ifstream& din, ifstream& ein);
void load_array_from_file(long double matrixelems[][N2], ifstream& fin);
void generate_lhn_matrix(long double omega, ofstream& fout);

#endif
