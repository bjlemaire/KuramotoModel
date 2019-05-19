#include "rk4_solver.hh"
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <omp.h>

int main(){
  double T_final = 1000.0;
  double hi_step = 1;
  double tolerance = 0.000001;
  double J = -0.5;
  double K = -0.5;
  double N_intrvls = 200;
  rk4 myRk4(hi_step, tolerance, J, K, N_intrvls);
  double t1 = omp_get_wtime();
  myRk4.compute_solution(T_final);
  double t2 = omp_get_wtime();
  myRk4.terminate();
  printf("Total Computation Time: %f\n", t2-t1);
  return 0;
}
