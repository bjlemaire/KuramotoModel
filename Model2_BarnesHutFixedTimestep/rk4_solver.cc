#include "rk4_solver.hh"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <omp.h>
#include <cfloat>

/** Constructor - variable initialization and Butcher tables. */
rk4::rk4(const int N_, double hi_step, double tol, double J_, 
        double K_, int n_intvls_, double barnes_theta_, bool enable_BH_)
        : enable_BH(enable_BH_), N(N_), AA(1.), BB(1.), J(J_), K(K_), 
        h_step(hi_step), last_h_step(hi_step),Atol(tol), Rtol(tol),
        fac(0.9), facmin(1./3.), facmax(3.), barnes_theta(barnes_theta_),
        SIZE_LIST(11), n_intvls(n_intvls_), not_found(1), step_counter(0), 
        A(new double[10]), B(new double[9]), C(new double[5]), 
        x0(new double[N]), y0(new double[N]), theta0(new double[N]), 
        vx0(new double[N]), vy0(new double[N]), omega0(new double[N]),
        x1(new double[N]), y1(new double[N]), theta1(new double[N]),
        x1h(new double[N]), y1h(new double[N]), theta1h(new double[N]),
        xlim(new double[3]), ylim(new double[3]), 
        xlim_next(new double[3]), ylim_next(new double[3]) {

          // Constructing the Butcher tables
          for(int i=0; i<10; i++){
            switch(i){
              case 0: A[i]=1./3.;  B[i]=1./8.; C[i]=0.;    break;
              case 1: A[i]=-1./3.; B[i]=3./8.; C[i]=1./3.; break;
              case 2: A[i]=1.;     B[i]=3./8.; C[i]=2./3.; break;
              case 3: A[i]=1.;     B[i]=1./8.; C[i]=1.;    break;
              case 4: A[i]=-1.;    B[i]=1./12.;C[i]=1.;    break;
              case 5: A[i]=1.;     B[i]=1./2.; break;
              case 6: A[i]=1./8.;  B[i]=1./4.; break;
              case 7: A[i]=3./8.;  B[i]=0.;    break;
              case 8: A[i]=3./8.;  B[i]=1./6.; break;
              case 9: A[i]=1./8.;  break;
              default: break;
            }
          }
        }

/** Destructor. */
rk4::~rk4(){
  zap(x0); zap(y0); zap(theta0); zap(x1); zap(y1); zap(theta1);
  zap(x1h); zap(y1h); zap(theta1h); zap(vx0); zap(vy0); zap(omega0);
  zap(A); zap(B); zap(C);
  zap(xlim); zap(ylim); zap(xlim_next); zap(ylim_next);
}

/** Given a particle and a root, find which square (or child node) 
  * the particle belongs to. */
void rk4::find_square(double x, double y, bool nlim){
  bool y_geq = y>ylim[1];
  bool x_leq = x<xlim[1];

  // Corresponding child node of the swarmalator - the location
  // value refers to the position of the child node in the node structure (see report
  // for more explanation).
  location = ( y_geq &&  x_leq)*7
            +( y_geq && !x_leq)*8
            +(!y_geq &&  x_leq)*9
            +(!y_geq && !x_leq)*10;

  // nlim: whether or not we want to update xlim and ylim (shrinking to 
  // the next subdomain).
  xlim_next[1] = nlim*( (  x_leq *0.5*(xlim[0]+xlim[1]))
                       +((!x_leq)*0.5*(xlim[1]+xlim[2])))
                + !nlim * xlim_next[1];
  xlim_next[0] = nlim*((!x_leq)?xlim[1]:xlim_next[0]) + !nlim * xlim_next[0];
  xlim_next[2] = nlim*( x_leq  ?xlim[1]:xlim_next[2]) + !nlim * xlim_next[2];
  ylim_next[1] = nlim*( (  y_geq *0.5*(ylim[1]+ylim[2]))
                       +((!y_geq)*0.5*(ylim[0]+ylim[1])))
                + !nlim * ylim_next[1];
  ylim_next[0] = nlim*(  y_geq ?ylim[1]:ylim_next[0]) + !nlim * ylim_next[0];
  ylim_next[2] = nlim*((!y_geq)?ylim[1]:ylim_next[2]) + !nlim * ylim_next[2];
}

/** Create a new node in the Barnes-Hut tree, given the data of a specific swarmalator. */
inline void rk4::push_node(int i, double x_, double y_, double sint_, double cost_){
  // Every new pushed node is a value node.
  barnes_list.push_back(1);

  // Decided that every value node has only one "children": itself.
  barnes_list.push_back(1);
  barnes_list.push_back(x_);
  barnes_list.push_back(y_);
  barnes_list.push_back(sint_);
  barnes_list.push_back(cost_);
  barnes_list.push_back(i);
  barnes_list.push_back(-1);
  barnes_list.push_back(-1);
  barnes_list.push_back(-1);
  barnes_list.push_back(-1);
}

/** Initialize the limits xlim and ylim as the smallest square
  * that includes all the swarmalators. */
inline void rk4::init_lims(){
  double low_  = fmin(minx-0.01*(maxx-minx), miny-0.01*(maxy-miny));
  double high_ = fmax(maxx+0.01*(maxx-minx), maxy+0.01*(maxy-miny));
  xlim[0] = low_;
  xlim[1] = 0.5*(low_+high_);
  xlim[2] = high_;
  ylim[0] = low_;
  ylim[1] = 0.5*(low_+high_);
  ylim[2] = high_;
  for (int i=0;i<3;i++){
    xlim_next[i] = xlim[i];
    ylim_next[i] = ylim[i];
  }
}

/** Perform the Barnes-Hut recursive algorithm on a given node - if the node
  * is a value, add its contribution to the net force. If the node is a root,
  * see report for algorithm.*/
void rk4::barnes_compute(int cidx_, int &i, double xi, double yi, double thi,
                         double &sumx, double &sumy, double &sumtheta, int &N_comp, 
                         double lgth){
  double mx=barnes_list[cidx_+2], my=barnes_list[cidx_+3],
         msin=barnes_list[cidx_+4], mcos=barnes_list[cidx_+5];
  int val = (int) barnes_list[cidx_], nchd=(int) barnes_list[cidx_+1],
      pid = (int) barnes_list[cidx_+6], n1=(int) barnes_list[cidx_+7],
      n2=(int) barnes_list[cidx_+8], n3=(int) barnes_list[cidx_+9],
      n4=(int) barnes_list[cidx_+10];
  double norm2 = (xi-mx)*(xi-mx) + (yi-my)*(yi-my);
  double norm = sqrt(norm2);

  // If node is a value node.
  if (val==1){
    if (i != pid){
      sumx += ((mx-xi)/norm)*(AA+J*(cos(thi)*mcos+sin(thi)*msin)) - BB *((mx-xi)/norm2);
      sumy += ((my-yi)/norm)*(AA+J*(cos(thi)*mcos+sin(thi)*msin)) - BB *((my-yi)/norm2);
      sumtheta += (msin*cos(thi)-sin(thi)*mcos)/norm;
      N_comp += nchd;
    }
  }
  else{

    // If d/s>Θ_bh, recursive algorithm on the 4 child nodes.
    if ((lgth/norm)>barnes_theta){
      if (n1!=-1) barnes_compute(n1,i,xi,yi,thi,sumx,sumy,sumtheta,N_comp, lgth/2.);
      if (n2!=-1) barnes_compute(n2,i,xi,yi,thi,sumx,sumy,sumtheta,N_comp, lgth/2.);
      if (n3!=-1) barnes_compute(n3,i,xi,yi,thi,sumx,sumy,sumtheta,N_comp, lgth/2.);
      if (n4!=-1) barnes_compute(n4,i,xi,yi,thi,sumx,sumy,sumtheta,N_comp, lgth/2.);
    }

    // If d/s<Θ_bh, add contribution of the node to the net force.
    else{
      sumx += nchd*(((mx-xi)/norm)*(AA+J*(cos(thi)*mcos+sin(thi)*msin)) - BB *((mx-xi)/norm2));
      sumy += nchd*(((my-yi)/norm)*(AA+J*(cos(thi)*mcos+sin(thi)*msin)) - BB *((my-yi)/norm2));
      sumtheta += nchd*(msin*cos(thi)-sin(thi)*mcos)/norm;
      N_comp += nchd;
    }
  }
}


/** Compute the N-body interactions using the Barnes-Hut algorithm. */
void rk4::smart_compute_xx(double t_, double* x_, double* y_, double* theta_, double* outputX, double* outputY, double* output_theta){

  // Create a fresh Barnes-Hut list
  barnes_list.clear();
  push_node(0, x_[0],y_[0],sin(theta_[0]),cos(theta_[0]));
  barnes_list[0] = 0;
  barnes_list[(y_[0]>0)?((x_[0]<0)?7:8):((x_[0]<0)?9:10)] = barnes_list.size();
  cidx = barnes_list.size();
  push_node(0, x_[0],y_[0],sin(theta_[0]),cos(theta_[0]));
  minx=DBL_MAX;
  maxx=-DBL_MAX;
  miny=DBL_MAX;
  maxy=-DBL_MAX;
  for (int i=0; i<N; i++){
    minx = (x_[i]<minx)?x_[i]:minx;
    maxx = (x_[i]>maxx)?x_[i]:maxx;
    miny = (y_[i]<miny)?y_[i]:miny;
    maxy = (y_[i]>maxy)?y_[i]:maxy;
  }
  // Initialize xlim, ylim, xlim_next, ylim_next
  init_lims();

  for(int i=1; i<N; i++){
    not_found = 1;
    cidx = 0;
    location = 0;
    init_lims();
    do{
      // If the current node is a root :
      if (barnes_list[cidx] == 0){

        // Find the child node corresponding to our particle i
        find_square(x_[i],y_[i],1);
        fidx = cidx+location;

        // If child node is empty node, create the new node.
        if (barnes_list[fidx]==-1){
          barnes_list[fidx]=barnes_list.size();
          push_node(i, x_[i],y_[i],sin(theta_[i]),cos(theta_[i]));

          // We increment the number of children of this node
          barnes_list[cidx+1] += 1;
          N_child = barnes_list[cidx+1];
          barnes_list[cidx+2] = ((N_child-1)*barnes_list[cidx+2] + x_[i])/N_child;
          barnes_list[cidx+3] = ((N_child-1)*barnes_list[cidx+3] + y_[i])/N_child;
          barnes_list[cidx+4] = ((N_child-1)*barnes_list[cidx+4] + sin(theta_[i]))/N_child;
          barnes_list[cidx+5] = ((N_child-1)*barnes_list[cidx+5] + cos(theta_[i]))/N_child;

          // Declare having found final node of particle i
          // - can leave the do-while loop and go to next
          // particle i+1.
          not_found = 0;
        }
        // If child node not empty, navigate to this next node.
        else {
          barnes_list[cidx+1] += 1;
          N_child = barnes_list[cidx+1];
          barnes_list[cidx+2] = ((N_child-1)*barnes_list[cidx+2] + x_[i])/N_child;
          barnes_list[cidx+3] = ((N_child-1)*barnes_list[cidx+3] + y_[i])/N_child;
          barnes_list[cidx+4] = ((N_child-1)*barnes_list[cidx+4] + sin(theta_[i]))/N_child;
          barnes_list[cidx+5] = ((N_child-1)*barnes_list[cidx+5] + cos(theta_[i]))/N_child;
          cidx = barnes_list[fidx];
          xlim[0]=xlim_next[0];
          xlim[1]=xlim_next[1];
          xlim[2]=xlim_next[2];
          ylim[0]=ylim_next[0];
          ylim[1]=ylim_next[1];
          ylim[2]=ylim_next[2];
        }
      }

      // If current node just a value : need to turn the node into a root.
      else{

        // Copy the values of value-node, and push them as new node
        ox = barnes_list[cidx+2];
        oy = barnes_list[cidx+3];
        osint = barnes_list[cidx+4];
        ocost = barnes_list[cidx+5];
        oid = barnes_list[cidx+6];
        barnes_list[cidx] = 0;
        find_square(ox,oy,0);
        barnes_list[cidx+location]=barnes_list.size();
        push_node(oid, ox, oy, osint, ocost);
        // Note: current index "cidx" remains identical.
      }
    } while(not_found);
  }

#pragma omp parallel for
  for(int i=0; i<N; i++){
    cidx = 0;
    outputX[i] = 0.;
    outputY[i] = 0.;
    output_theta[i] = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumtheta = 0.;
    int N_comp = 0;
    barnes_compute(0, i, x_[i], y_[i], theta_[i], sumx, sumy, sumtheta, N_comp, maxx-minx);
    outputX[i] = (1./float(N_comp))*sumx;
    outputY[i] = (1./float(N_comp))*sumy;
    output_theta[i] = (float(K)/float(N_comp))*sumtheta;
  }

  barnes_list.clear();
}

/** Compute the N-body interactions using a "brute-force" calculation. */
void rk4::compute_xx(double t_, double* x_, double* y_, double* theta_, double* outputX, double* outputY, double* output_theta){

  // Calculation of the net force exerted on each swarmalator.
#pragma omp parallel for
  for (int i=0; i<N; i++){

    // Array initialization
    outputX[i] = 0.;
    outputY[i] = 0.;
    output_theta[i] = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumtheta = 0.;

    // Brute force calculation of the N-body interactions
    for (int j=0; j<N; j++){
      if (j!=i){
        double norm2 = (x_[j]-x_[i])*(x_[j]-x_[i]) + (y_[j]-y_[i])*(y_[j]-y_[i]);
        double norm = sqrt(norm2);
        sumx += ((x_[j]-x_[i])/norm)*(AA+J*cos(theta_[j]-theta_[i])) - BB*((x_[j]-x_[i])/norm2);
        sumy += ((y_[j]-y_[i])/norm)*(AA+J*cos(theta_[j]-theta_[i])) - BB*((y_[j]-y_[i])/norm2);
        sumtheta += sin(theta_[j]-theta_[i])/norm;
      }
    }
    outputX[i] += (1./float(N))*sumx;
    outputY[i] += (1./float(N))*sumy;
    output_theta[i] += (float(K)/float(N))*sumtheta;
  }
}

/** Compute the G vectors (see FSAL RK4 method, AM225 lecture slides). */
void rk4::compute_Gs(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta){
  // First, we calculate G1:
#pragma omp parallel for
  for(int i=0; i<(N); i++){
    Gs_x[i]=x0[i];
    Gs_x[i+1*N]=x0[i];
    Gs_x[i+2*N]=x0[i];
    Gs_x[i+3*N]=x0[i];
    Gs_x[i+4*N]=x0[i];
    Gs_y[i]=y0[i];
    Gs_y[i+1*N]=y0[i];
    Gs_y[i+2*N]=y0[i];
    Gs_y[i+3*N]=y0[i];
    Gs_y[i+4*N]=y0[i];
    Gs_theta[i]=theta0[i];
    Gs_theta[i+1*N]=theta0[i];
    Gs_theta[i+2*N]=theta0[i];
    Gs_theta[i+3*N]=theta0[i];
    Gs_theta[i+4*N]=theta0[i];
  }
  // Then, we compute f(G1):
  if (enable_BH) smart_compute_xx(t+C[0]*h_step, x0, y0, theta0, ff_x, ff_y, ff_theta);
  else compute_xx(t+C[0]*h_step, x0, y0, theta0, ff_x, ff_y, ff_theta);

  // Calculating G2:
#pragma omp parallel for
  for(int i=0; i<N; i++){
    Gs_x[i+1*N]     += A[0]*h_step*ff_x[i];
    Gs_y[i+1*N]     += A[0]*h_step*ff_y[i];
    Gs_theta[i+1*N] += A[0]*h_step*ff_theta[i];
  }

  // Computing f(G2):
  if (enable_BH) smart_compute_xx(t+C[1]*h_step, (Gs_x+1*N), (Gs_y+1*N), (Gs_theta+1*N), (ff_x+1*N), (ff_y+1*N), (ff_theta+1*N));
  else compute_xx(t+C[1]*h_step, (Gs_x+1*N), (Gs_y+1*N), (Gs_theta+1*N), (ff_x+1*N), (ff_y+1*N), (ff_theta+1*N));

  // Calculating G3:
  for (int i=0; i<N; i++){
    Gs_x[i+2*N]     += A[1]*h_step*ff_x[i];
    Gs_x[i+2*N]     += A[2]*h_step*ff_x[i+1*N];
    Gs_y[i+2*N]     += A[1]*h_step*ff_y[i];
    Gs_y[i+2*N]     += A[2]*h_step*ff_y[i+1*N];
    Gs_theta[i+2*N] += A[1]*h_step*ff_theta[i];
    Gs_theta[i+2*N] += A[2]*h_step*ff_theta[i+1*N];
  }

  // Computing f(G3):
  if (enable_BH) smart_compute_xx(t+C[2]*h_step, (Gs_x+2*N), (Gs_y+2*N), (Gs_theta+2*N), (ff_x+2*N), (ff_y+2*N), (ff_theta+2*N));
  else compute_xx(t+C[2]*h_step, (Gs_x+2*N), (Gs_y+2*N), (Gs_theta+2*N), (ff_x+2*N), (ff_y+2*N), (ff_theta+2*N));


  // Calculating G4:
#pragma omp parallel for
  for (int i=0; i<N; i++){
    Gs_x[i+3*N]     += A[3]*h_step*ff_x[i];
    Gs_x[i+3*N]     += A[4]*h_step*ff_x[i+1*N];
    Gs_x[i+3*N]     += A[5]*h_step*ff_x[i+2*N];
    Gs_y[i+3*N]     += A[3]*h_step*ff_y[i];
    Gs_y[i+3*N]     += A[4]*h_step*ff_y[i+1*N];
    Gs_y[i+3*N]     += A[5]*h_step*ff_y[i+2*N];
    Gs_theta[i+3*N] += A[3]*h_step*ff_theta[i];
    Gs_theta[i+3*N] += A[4]*h_step*ff_theta[i+1*N];
    Gs_theta[i+3*N] += A[5]*h_step*ff_theta[i+2*N];
  }

  // Computing f(G4):
  if (enable_BH) smart_compute_xx(t+C[3]*h_step, (Gs_x+3*N), (Gs_y+3*N), (Gs_theta+3*N), (ff_x+3*N), (ff_y+3*N), (ff_theta+3*N));
  else compute_xx(t+C[3]*h_step, (Gs_x+3*N), (Gs_y+3*N), (Gs_theta+3*N), (ff_x+3*N), (ff_y+3*N), (ff_theta+3*N));


  // Calculating G5:
#pragma omp parallel for
  for (int i=0; i<N; i++){
    Gs_x[i+4*N]     += A[6]*h_step*ff_x[i];
    Gs_x[i+4*N]     += A[7]*h_step*ff_x[i+1*N];
    Gs_x[i+4*N]     += A[8]*h_step*ff_x[i+2*N];
    Gs_x[i+4*N]     += A[9]*h_step*ff_x[i+3*N];
    Gs_y[i+4*N]     += A[6]*h_step*ff_y[i];
    Gs_y[i+4*N]     += A[7]*h_step*ff_y[i+1*N];
    Gs_y[i+4*N]     += A[8]*h_step*ff_y[i+2*N];
    Gs_y[i+4*N]     += A[9]*h_step*ff_y[i+3*N];
    Gs_theta[i+4*N] += A[6]*h_step*ff_theta[i];
    Gs_theta[i+4*N] += A[7]*h_step*ff_theta[i+1*N];
    Gs_theta[i+4*N] += A[8]*h_step*ff_theta[i+2*N];
    Gs_theta[i+4*N] += A[9]*h_step*ff_theta[i+3*N];
  }

  // Computing f(G5):
  if (enable_BH) smart_compute_xx(t+C[4]*h_step, (Gs_x+4*N), (Gs_y+4*N), (Gs_theta+4*N), (ff_x+4*N), (ff_y+4*N), (ff_theta+4*N));
  else compute_xx(t+C[4]*h_step, (Gs_x+4*N), (Gs_y+4*N), (Gs_theta+4*N), (ff_x+4*N), (ff_y+4*N), (ff_theta+4*N));
}

/** Compute y1 and y1h using the Butcher tables and 
  * the FSAL RK4 numerical integration scheme. */
void rk4::compute_y1y1h(double t, double* Gs_x, double* ff_x, double* Gs_y, double* ff_y, double* Gs_theta, double* ff_theta){

  // Initialize the Scaling arrays.
  double* sc_x = new double[N];
  double* sc_y = new double[N];
  double* sc_theta = new double[N];

#pragma omp parallel for
  for (int i=0; i<N; i++){
    sc_x[i] = 0.;
    sc_y[i] = 0.;
    sc_theta[i] = 0.;
    x1[i] = 0.;
    y1[i] = 0.;
    theta1[i] = 0.;
    x1h[i] = 0.;
    y1h[i] = 0.;
    theta1[i] = 0.;
  }

  // Calculate x1, y1, Θ1, and their "hat" values (see FSAL RK4 integration scheme).
#pragma omp parallel for
  for (int i=0; i<N; i++){
    x1[i]      = x0[i];
    y1[i]      = y0[i];
    theta1[i]  = theta0[i];
    x1h[i]     = x0[i];
    y1h[i]     = y0[i];
    theta1h[i] = theta0[i];
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<4; j++){
      x1[i] += B[j]*h_step*ff_x[i+j*N];
      y1[i] += B[j]*h_step*ff_y[i+j*N];
      theta1[i] += B[j]*h_step*ff_theta[i+j*N];
    }
  }

  // Deallocating arrays.
  zap(sc_x); zap(sc_y); zap(sc_theta);
}

/** Performs 3rd-order Hermite polynomial interpolation to provide dense output. */
void rk4::hermite(double actual_t, double myTheta, char* filenameDense){
  std::ofstream myDense;
  myDense.open(filenameDense);
  double* f0_x = new double[N];
  double* f0_y = new double[N];
  double* f0_theta = new double[N];
  double* f1_x = new double[N];
  double* f1_y = new double[N];
  double* f1_theta = new double[N];

  // If enable_BH is true, perform integration using Barnes-Hut algorithm.
  if (enable_BH){
    smart_compute_xx(actual_t, x0, y0, theta0, f0_x, f0_y, f0_theta);
    smart_compute_xx(actual_t + last_h_step, x1, y1, theta1, f1_x, f1_y, f1_theta);
  }

  // Else, perform integration using brute-force calculation.
  else {
    compute_xx(actual_t, x0, y0, theta0, f0_x, f0_y, f0_theta);
    compute_xx(actual_t + last_h_step, x1, y1, theta1, f1_x, f1_y, f1_theta);
  }

  // Provide dense output.
  for(int i=0; i<N; i++){
    double u_x = (1-myTheta)*x0[i] + myTheta*x1[i] + myTheta*(myTheta-1)*((1-2*myTheta)*(x1[i]-x0[i]) + (myTheta-1)*last_h_step*f0_x[i]+myTheta*last_h_step*f1_x[i]);
    double u_y = (1-myTheta)*y0[i] + myTheta*y1[i] + myTheta*(myTheta-1)*((1-2*myTheta)*(y1[i]-y0[i]) + (myTheta-1)*last_h_step*f0_y[i]+myTheta*last_h_step*f1_y[i]);
    double u_theta = (1-myTheta)*theta0[i] + myTheta*theta1[i] + myTheta*(myTheta-1)*((1-2*myTheta)*(theta1[i]-theta0[i]) + (myTheta-1)*last_h_step*f0_theta[i]+myTheta*last_h_step*f1_theta[i]);
    myDense<<(actual_t + myTheta*last_h_step)<<" "<<u_x <<" "<< u_y<<" "<<u_theta<<std::endl;
  }
  myDense.close();

  // Array deallocation.
  zap(f0_x); zap(f0_y); zap(f0_theta);
  zap(f1_x); zap(f1_y); zap(f1_theta);
}

/** Provide a dense output. */
void rk4::dense_output(double t_){
  int m1=0;
  if (float(int(t_/dense_stpsze)) == t_/dense_stpsze){
    m1 = int(t_/dense_stpsze)-1;
  }
  else{
    m1 = int(t_/dense_stpsze);
  }
  int m2 = int((t_+last_h_step)/dense_stpsze);
  for (int k=(m1+1); k<(m2+1); k++){
    char filenameDense[32];
    sprintf(filenameDense, "Dense%04d.txt", k);
    hermite(t_, (dense_stpsze*k-t_)/last_h_step, filenameDense);
  }
}

/** Prepare all the variables to pass their values to reinitialize the system. */
void rk4::nextStep(){
  for(int i=0; i<N; i++){
    x0[i]=x1[i];
    y0[i]=y1[i];
    theta0[i]=theta1[i];
  }
}

/** Perform the calculation from t=0 to t=T_final. */
void rk4::compute_solution(double T_final_){
  T_final = T_final_;
  dense_stpsze = double(T_final/double(n_intvls));
  double t=0;
  double* Gs_x = new double[5*N];
  double* ff_x = new double[5*N];
  double* Gs_y = new double[5*N];
  double* ff_y = new double[5*N];
  double* Gs_theta = new double[5*N];
  double* ff_theta = new double[5*N];
  initialize();
  while(t<T_final){
    compute_Gs(t, Gs_x, ff_x, Gs_y, ff_y, Gs_theta, ff_theta);
    compute_y1y1h(t, Gs_x, ff_x, Gs_y, ff_y, Gs_theta, ff_theta);
    dense_output(t);
    t += last_h_step;
    //printf("t=%f, next_step=%f\n", t, h_step);
    step_counter += 1;
    nextStep();
  };
  //printf("Done! It tool us %d steps to perform the entire integration.", step_counter);
  zap(Gs_x); zap(ff_x); zap(Gs_y); zap(ff_y);
  zap(Gs_theta); zap(ff_theta);
}

/** Uniformly distribute x,y within a disk of radius 1 centered at the origin
  * Uniformly distribute Θ between -π and +π
  * Set the values of Vx, Vy, and ω to the ones defined by the author of the code. */
void rk4::initialize(){
  srand (static_cast <unsigned> (time(0)));
  float max_xy = 2;
  float max_angle = 2*M_PI;
  for(int i=0; i<N; i++){
    float x;
    float y;
    float phse;
    do {
      x = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/max_xy));
      y = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/max_xy));
      phse = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/max_angle));
      x-=1;
      y-=1;
    } while ((x*x + y*y)>1);
    x0[i] = x;
    y0[i] = y;
    theta0[i] = phse;
    vx0[i] = 0;
    vy0[i] = 0;
    omega0[i] = 0;
  }
}

