#ifndef _FUNCTIONS_H_INCLUDED_
#define _FUNCTIONS_H_INCLUDED_

#include <iostream>
#include <cmath>
#include <chrono>
#include "common_globals.h"

using namespace std;
using namespace std::chrono;

const int N5=10;

double psi(int i, double y);
double integrate_y(double (*func)(double, int, int, int, int), int l, int p, int n, int q ); //func = v_d/_ex
double integrate_x(double y, double (*func_x)(double, double,int, int), int l, int n); //func_x= v_x
double v_x(double x, double y, int a, int b);
double v_d(double y, int l, int p, int n, int q);
double v_ex(double y, int l, int p, int n, int q);
string current_time_str(void);
void createfilename(string& filename, string option, int number_of_mesh, double omega);
void load_array_from_file(double direct[][N5][N5][N5],double exchange[][N5][N5][N5], ifstream& din, ifstream& ein);
void load_array_from_file(double matrixelems[][N5], ifstream& fin);
void generate_lhn_matrix(double omega, ofstream& fout);
double integrate_x2(double (*func)(double, int, int), int l, int n);
void show_time(milliseconds begin_ms, milliseconds end_ms);


#endif
