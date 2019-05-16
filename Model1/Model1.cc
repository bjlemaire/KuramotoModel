#include "rk4_solver.hh"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>

int main(){
  double T_final = 100.0;
  double hi_step = 1;
  double tolerance = 0.000001;
  double J = 0.1;
  double K = 1;
  double N_intrvls = 2;
  rk4 myRk4(hi_step, tolerance, J, K, N_intrvls);
  myRk4.compute_solution(T_final);
  myRk4.terminate();
  return 0;
}
