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
  rk4 myRk4(hi_step, tolerance, J, K, N_intrvls,0.5);
  double xxx[4] = {-0.8, 0.8, 0.9, 0.55};
  double yyy[4] = {-0.9, 0.8, 0.95, 0.65};
  double ttt[4] = {1.5,2.6,3.5,4.5};
  double optx[4];
  double opty[4];
  double optth[4];
  myRk4.smart_compute_xx(3,xxx,yyy,ttt,optx,opty,optth);
  //double t1 = omp_get_wtime();
  //myRk4.compute_solution(T_final);
  //double t2 = omp_get_wtime();
  //myRk4.terminate();
  //printf("Total Computation Time: %f\n", t2-t1);
  return 0;
}
