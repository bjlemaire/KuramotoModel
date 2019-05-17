#ifndef RK4_SOLVER_HH
#define RK4_SOLVER_HH

#include <iostream>
#include <cassert>
#include <vector>
//#define N 1250

class rk4 {
  public:

    // Bool: if 0, use the initial code implementation
    //       if 1, use the Barnes-Hut code implementation
    bool enable_BH;

    // AA, BB: rescaling parameters from the Kuramoto model equations.
    // K, J: set of parameters of the Kuramoto model.
    // h_step, last_h_step: timestep sizes.
    // Atol, Rtol: tolerance thresholds of the adaptive timestep method.
    // fac, facmin, facmax: parameters of the adaptive timestep method.
    // T_final: final simulation time.
    // dense_stpsze: timestep size necessary to obtain the desired dense output
    // (obtained from T_final and the number of interval input from the operator)
    // barnes_theta: the theta threshold (Barnes-Hut Threshold) for the BH algorithm.
    // minx, maxx, miny, maxy: spatial limits of the N swarmalators.
    // ox, oy, osint, ocost, oid: variables used to turn a "Value" node into 
    // a "Root" node and a "Value" node.
    double AA, BB, J, K, h_step, last_h_step, Atol, Rtol, fac, facmin, facmax, T_final,
           dense_stpsze, barnes_theta, minx, maxx, miny, maxy, ox, oy, osint, ocost, oid;
    // N: number of swarmalators.
    // n_intvls: number of intervals for the dense output.
    // not_found: 
    // cidx: Current Index
    // fidx: Futur Index
    int N, n_intvls, not_found, location, cidx, fidx, N_child, step_counter, SIZE_LIST;
    double *A, *B, *C, *sc, *x0, *y0, *theta0, *vx0, *vy0, *omega0, *x1, *y1, *theta1,
           *x1h, *y1h, *theta1h, *vx1, *vy1, *omega1, *sc_x, *sc_y, *sc_theta, *bnodes_idx,
           *xlim, *ylim, *xlim_next, *ylim_next;
    std::vector<double> barnes_list;
    std::vector<double> bnodes;
    rk4(int N_, double hi_step, double tol, double J_, double K_, int n_intvls_, double barnes_theta_, bool enable_BH_);
    ~rk4();
    void initialize();
    void compute_solution(double T_final_);
  private:
    void compute_xx(double t_, double* x_, double* y_, double* theta_, double* outputX, double* outputY, double* output_theta);
    void compute_Gs(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta);
    void compute_y1y1h(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta);
    void dense_output(double t_);
    void hermite(double actual_t, double myTheta, char* filenameDense);
    void nextStep();
    void zap(double* myArray){
      {assert(myArray!=NULL);}
      delete [] myArray;
      myArray = NULL;
    };
    void barnes_compute(int cidx_, int &i, double xi, double yi, double thi,
                         double &sumx, double &sumy, double &sumtheta, int &N_comp, double lgth);
    void smart_compute_xx(double t_, double* x_, double* y_, double* theta_,
                          double* outputX, double* outputY, double* output_theta );
    inline void init_lims();
    void find_square(double x, double y, bool nlim);
    inline void push_node(int i, double x_, double y_, double sint_, double cost_);


};

#endif
