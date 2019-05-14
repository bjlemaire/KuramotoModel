#include "rk4_solver.hh"
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <omp.h>

int main(){
  double T_final = 1.0;
  double hi_step = 1;
  double tolerance = 0.000001;
  double J = 1;
  double K = 1;
  double N_intrvls = 2;
  rk4 myRk4(hi_step, tolerance, J, K, N_intrvls,0.9);
  //double xxx[3] = {-0.8, 0.8, 0.9};
  //double yyy[3] = {-0.9, 0.8, 0.95};
  //double ttt[3] = {1.5,2.6,3.5};
  //double optx[3];
  //double opty[3];
  //double optth[3];
  //myRk4.smart_compute_xx(3,xxx,yyy,ttt,optx,opty,optth);
  double t1 = omp_get_wtime();
  myRk4.compute_solution(T_final);
  double t2 = omp_get_wtime();
  myRk4.terminate();
  printf("Total Computation Time: %f\n", t2-t1);
  return 0;
}
