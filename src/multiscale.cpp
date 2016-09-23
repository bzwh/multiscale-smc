#include <iostream>
#include <ctime>
#include <sstream>
#include <gsl/gsl_rng.h>
#include "Farm.hpp"
#include "Smc.hpp"

#include <omp.h>

using namespace std;


int main ( int argc, char* argv[])  {
  const int nthrds = 48;
  const int nparts = 5000;
  const int nround = 20;


  Farms frms;     // Initialise and  load outside of parallel region.
  frms.bigload(); // Nothing inside changes, only need one copy - slows down when threads race
  Smc smc(nthrds,nparts,nround,frms);

  smc.run();


  return(0);
}












