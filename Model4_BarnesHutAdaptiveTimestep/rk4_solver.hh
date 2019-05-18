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

    // Total number of swarmalators
    const int N;

    // -  AA, BB:   rescaling parameters from the Kuramoto model equations.
    // -  K, J:     set of parameters of the Kuramoto model.
    // -  h_step,   last_h_step:  timestep sizes.
    // -  Atol, Rtol:           tolerance thresholds of the adaptive timestep method.
    // -  fac, facmin, facmax:  parameters of the adaptive timestep method.
    // -  T_final:  final simulation time.
    // -  dense_stpsze: timestep size necessary to obtain the desired dense output
    //                  (obtained from T_final and the number of interval input 
    //                  from the operator).
    // -  barnes_theta: theta threshold (Barnes-Hut Threshold) for the BH algorithm.
    // -  minx, maxx, miny, maxy:     spatial limits of the N swarmalators.
    // -  ox, oy, osint, ocost, oid:  variables used to turn a "Value" node into 
    //                                a "Root" node and a "Value" node.
    double AA, BB, J, K, h_step, last_h_step, Atol, Rtol, fac, facmin, facmax, T_final,
           dense_stpsze, barnes_theta, minx, maxx, miny, maxy, ox, oy, osint, ocost, oid;


    // Size of the data structure of 1 node of the Barnes-Hut Tree
    const int SIZE_LIST;

    // N: number of swarmalators.
    // n_intvls: number of intervals for the dense output.
    // not_found: 
    // cidx: Current Index
    // fidx: Futur Index
    int n_intvls, not_found, location, cidx, fidx, N_child, step_counter;
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
    // Compute the interaction in the system of 
    // equation (right terms) using brute force calculation.
    void compute_xx(double t_, double* x_, double* y_, double* theta_, double* outputX, double* outputY, double* output_theta);

    // Compute the G vectors in the FSAL RK4 scheme.
    void compute_Gs(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta);

    // Compute y1 and y1h using the Butcher tables (RK4 scheme)
    void compute_y1y1h(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta);

    // Provide a dense output
    void dense_output(double t_);

    // Perform 3rd-order Hermite polynomial interpolation
    // to provide dense output.
    void hermite(double actual_t, double myTheta, char* filenameDense);

    // Prepare all the variables to pass their values to reinitialize the system.
    void nextStep();

    // Clean deallocation of the pointers.
    void zap(double* myArray){
      {assert(myArray!=NULL);}
      delete [] myArray;
      myArray = NULL;
    };

    // Recursively compute the interaction on a given node of the BH-Tree
    void barnes_compute(int cidx_, int &i, double xi, double yi, double thi,
                         double &sumx, double &sumy, double &sumtheta, int &N_comp, double lgth);

    // Use the Barnes-Hut algorithm to efficienttly compute of the N-body interactions.
    void smart_compute_xx(double t_, double* x_, double* y_, double* theta_,
                          double* outputX, double* outputY, double* output_theta );

    // Initialize xlim and ylim to the values of the first domain (smallest
    // square containing all the swarmalators).
    inline void init_lims();

    void find_square(double x, double y, bool nlim);
    inline void push_node(int i, double x_, double y_, double sint_, double cost_);


};

#endif
