#ifndef RK4_SOLVER_HH
#define RK4_SOLVER_HH

#include <iostream>
#include <cassert>
#include <vector>

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
    // -  barnes_theta: Θ threshold (Barnes-Hut Threshold) for the BH algorithm.
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

    // - A,B,C:             Butcher Tables used for FSAL Runge-Kutta 4th order numerical scheme.
    // - sc:                Used for adaptive method. See AM225 lecture slides.
    // - x0,y0,theta0:      Initial values for x,y,Θ (initial meaning before each new timestep).
    // - vx0, vy0, omega0:  Initial values for velocity along x, y, and ω.
    // - x1, y1, theta1:    Final values for x,y,Θ (final meaning after each new timestep).
    // - x1h, y1h, theta1h: Values for x1-hat, y1-hat, and Θ-hat used for adaptive timestep method.
    // - vx1, vy1, omega1:  Final values for velocity along x, y, and ω.
    //
    // - xlim, ylim, xlim_next, ylim_next: Smallest squares including all the swarmalators.
    //                                     Values along x and y before and after each timestep 
    //                                     respectively.
    double *A, *B, *C, *sc, *x0, *y0, *theta0, *vx0, *vy0, *omega0, *x1, *y1, *theta1,
           *x1h, *y1h, *theta1h, *vx1, *vy1, *omega1, *xlim, *ylim, *xlim_next, *ylim_next;

    // The Barnes-Hut Tree vector.
    std::vector<double> barnes_list;

    // Constructor.
    rk4(int N_, double hi_step, double tol, double J_, double K_, int n_intvls_, double barnes_theta_, bool enable_BH_);

    // Destructor.
    ~rk4();

    // Uniformly distribute x,y within a disk of radius 1 centered at the origin.
    // Uniformly distribute Θ between -π and π.
    // Set the values of Vx, Vy, ω to the ones defined by the author of the code.
    void initialize();

    // Perform the calculation from t=0 to t=T_final.
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

    // Given a particle and a Root node, find each square (or Child node) 
    // the particle belongs to.
    void find_square(double x, double y, bool nlim);

    // Creates a new Value node in the Barnes-Hut tree with the data passed as arguments.
    inline void push_node(int i, double x_, double y_, double sint_, double cost_);


};

#endif
