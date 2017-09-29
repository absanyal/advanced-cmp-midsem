#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>

#include "common-globals.h"
#include "functions.h"

using namespace std;

ifstream din("direct.txt");
ifstream ein("exchange.txt");

int main()
{
  long double direct[3][3][3][3];
  long double exchange[3][3][3][3];

  load_array_from_file(direct,exchange,din,ein);
}
