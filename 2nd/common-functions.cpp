#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "common_globals.h"

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace boost::math;


void show_time(milliseconds begin_ms, milliseconds end_ms)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout << "Diagonalization took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << "Diagonalization took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << "Diagonalization took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << "Diagonalization took " << t << "time. Wrong t received.\n"; }
}
