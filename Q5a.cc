#include "rk4_solver.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>

int main(){
  double T_final = 200.0;
  double hi_step = 1;
  double tolerance = 0.000001;
  double J = 1.;
  double K = -0.2;
  double N_intrvls = 600;
  rk4 myRk4(hi_step, tolerance, J, K, N_intrvls);
  std::cout<<"Ready to dive in!\n";
  myRk4.compute_solution(T_final);
  myRk4.terminate();
  return 0;
}