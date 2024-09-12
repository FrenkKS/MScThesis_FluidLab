#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "fluidlab_pack.h"

// Comment or uncomment the line below to do Oscillatory shear or start-up shear flow
#define OSCILLATORY 1

// Comment or uncomment these to use the analytically calculated second invariant and its derivative
// or the ones determined using finite differences. If ANALYTICAL_SECOND_INV_DERIV is True,
// ANALYTICAL should also be True
#define ANALYTICAL 1
#define ANALYTICAL_SECOND_INV_DERIV 1

// Function to determine the sign of the argument
#define sign_func(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))

// Define model parameters
double Re = 0.1; // Reynolds number
double regularization = 1e-3; // Regularization parameter

double gamma_0 = 5.29; // Imposed shear rate (amplitude if oscillatory shear)

double Wi = 0.3; // Weissenberg number
double Bi = 7.272; // Bingham number
double beta = 1.7; // viscosity ratio
double n = 0.416; // power index

int decs = 3; // Number of decimals to round Wi, Bi, and beta on

// Fields to store log-conform variables
scalar lambdav[], eta_p[];

// Field for storing second invariant values
scalar second_inv_field[];

// Field for storing positional derivatives of second invariant
scalar second_inv_deriv_x_field[], second_inv_deriv_y_field[];

// The top and bottom boundary conditions are those of a Couette flow.
#ifdef OSCILLATORY
  u.t[top] = dirichlet (y*cos(t/gamma_0));
  u.t[bottom] = dirichlet (y*cos(t/gamma_0));
#else
  u.t[top] = dirichlet (y);
  u.t[bottom] = dirichlet (y);
#endif

// Time to end the simulation
double final_time;

// Custom time variable for shear calculation
double time_variable;

int main (int argc, char * argv[]) {

  // Change dimensionless variables to appropriate values based on imposed shear rate
  Wi = Wi * gamma_0;
  Bi = Bi / gamma_0;
  beta = beta * pow(gamma_0, n - 1.0);

  double scale = pow(10, decs);
  Wi = round(Wi * scale) / scale;
  Bi = round(Bi * scale) / scale;
  beta = round(beta * scale) / scale;

  // Creating the domain (rectangle [-10, 10] x [-1, 1])
  // The extra y-coordiantes will be masked in the init
  size(20.0);
  origin(-10.0, -10.0);
  periodic(right);
  init_grid(256);

  // Because we access the momentum equation via two-phase.h,
  // we need to set both fluids to be the same
  rho1 = rho2 = Re;
  mu1 = mu2 = Wi;

  
  // The viscoelastic fields will be set below.
  mup = eta_p;
  lambda = lambdav;

  // Define starting time
  time_variable = 0.0;

  // Making the poisson solver tolerance a little bit tighter
  TOLERANCE = 1e-4;

  // Open simulation folders for storing the results
  // and define the simulation time
  double analytical_indicator;
  #ifdef ANALYTICAL
    analytical_indicator = 1.0;
  #else
    analytical_indicator = 0.0;
  #endif

  #ifdef OSCILLATORY
    OpenSimulationFolder("oscillatory_reg%g_Wi%g_Bi%g_beta%g_n%g_Re%g_amp%g_an%g", regularization, Wi, Bi, beta, n, Re, gamma_0, analytical_indicator);
    double number_periods = 1.0;
    final_time = number_periods*2.0*M_PI;
  #else
    OpenSimulationFolder("startup_reg%g_Wi%g_Bi%g_beta%g_n%g_Re%g_amp%g_an%g", regularization,  Wi, Bi, beta, n, Re, gamma_0, analytical_indicator);
    final_time = 10.0;
  #endif

  // Running the simulation
  run();

  // Closing the folder and open files
  CloseSimulationFolder();
}

// Set the maximum timestep
event control_timestep(i++) {
  DT = 0.001;
}

// Store the stress values for plotting purposes
scalar txx_old[], txy_old[], tyy_old[];
event logfile (i += 1; t<=final_time) {

  // Get the stress at the center of the domain
  double txx = interpolate(tau_p.x.x, 0.0, 0.0);
  double txy = interpolate(tau_p.x.y, 0.0, 0.0);
  double tyy = interpolate(tau_p.y.y, 0.0, 0.0);

  // Check how much the stress has changed since the last iteration
  // This can be used to notice diverging stresses early on
  double change_txx = change(tau_p.x.x, txx_old);
  double change_txy = change(tau_p.x.y, txy_old);
  double change_tyy = change(tau_p.y.y, tyy_old);
  double difference = change_txx + change_txy + change_tyy;

  // Printing the stresses to a log file
  printf("Time step (%d, %lf, %g): %e ... %e %e %e\n", i, t, dt, difference, txx, txy, tyy);
  PrintLog("%d %lf %e %e %e %e\n", i, t, difference, txx, txy, tyy);
}

event print_solution(t+=10.0) {
  // Printing the properties below into a VTK file. Can be viewed in Paraview
  scalar *list = {u.x, p, tau_p.x.x, tau_p.x.y, tau_p.y.y};
  const char *list_names[] = {"vel-u", "pressure", "taup_xx", "taup_xy", "taup_yy"};
  PrintMeshVTK_Binary_Float(i, t, list, list_names);
}

event init (t = 0) {

  // Mask the regios we want to exclude
  mask(y>=1.0 ? top : none);
  mask(y<=-1.0 ? bottom : none);

  // Initializing volume fractions (again because we are using two-phase)
  foreach() {
    f[] = ( y<1.0 && y>-1.0 ) ? 1.0 : 0.0;
    u.x[] = y;
  }

  // Set the initial value for the second derivative at each cell
  // This is necessary for the finite differences method for calculating its derivative
  foreach() {

    double second_inv;

    #ifdef ANALYTICAL

      double shear_rate;

      #ifdef OSCILLATORY
        shear_rate = cos(0);
      #else
        shear_rate = 1.0;
      #endif

      // Calculate second invariant for shear experiments analytically
      second_inv = sqrt(2) * fabs(shear_rate);
    
    #else
      double dux_dx = (u.x[] - u.x[-1,0]) / Delta;
      double dux_dy = (u.x[] - u.x[0,-1]) / Delta;
      double duy_dx = (u.y[] - u.y[-1,0]) / Delta;
      double duy_dy = (u.y[] - u.y[0,-1]) / Delta;

      // Calculate second invariant for 2D in general
      second_inv = sqrt(4 * pow(dux_dx, 2) + 4 * pow(duy_dy, 2) + 2 * pow(dux_dy + duy_dx, 2));
    #endif

    second_inv_field[] = second_inv;
  }
}

// Updating the EVP properties
event properties (i++) {

  // This if-statement is necessary to ensure correct calculations
  if (i > 0) {

    // Set the positional derivatives of the second invariant
    // This is necessary for the finite differences method for calculating its derivative
    foreach() {
      second_inv_deriv_x_field[] = (second_inv_field[] - second_inv_field[-1, 0]) / Delta;
      second_inv_deriv_y_field[] = (second_inv_field[] - second_inv_field[0, -1]) / Delta;
    }

    foreach() {

      double shear_rate;
      double shear_rate_rate;

      double second_inv;

      #ifdef ANALYTICAL

        #ifdef OSCILLATORY
          shear_rate = cos(time_variable / gamma_0);
          shear_rate_rate = - (1 / gamma_0) * sin(time_variable / gamma_0);
        #else
          shear_rate = 1.0;
          shear_rate_rate = 0.0;
        #endif

        // Calculate second invariant for shear experiments analytically
        second_inv = sqrt(2) * fabs(shear_rate);
      
      #else
        
        double dux_dx = (u.x[] - u.x[-1,0]) / Delta;
        double dux_dy = (u.x[] - u.x[0,-1]) / Delta;
        double duy_dx = (u.y[] - u.y[-1,0]) / Delta;
        double duy_dy = (u.y[] - u.y[0,-1]) / Delta;

        // Calculate second invariant for 2D in general
        second_inv = sqrt(4 * pow(dux_dx, 2) + 4 * pow(duy_dy, 2) + 2 * pow(dux_dy + duy_dx, 2));
      #endif

      double second_inv_deriv;

      #ifdef ANALYTICAL_SECOND_INV_DERIV
        // Calculate material derivative of second invariant for shear experiments analytically
        second_inv_deriv = sqrt(2) * sign_func(shear_rate) * shear_rate_rate;

      #else
        // Calculate material derivative of second invariant for 2D in general
        double second_inv_deriv_t = (second_inv - second_inv_field[]) / dt;
        second_inv_deriv = second_inv_deriv_t + u.x[] * second_inv_deriv_x_field[] + u.y[] * second_inv_deriv_y_field[];

      #endif

      second_inv_field[] = second_inv;

      double yield_term = Bi / (second_inv + regularization);

      double HB_term;

      if(n < 1.0) {
        HB_term = (beta * pow(second_inv, n)) / (second_inv + regularization);
      }
      else {
        HB_term = beta * pow(second_inv, n - 1.0);
      }

      double s = yield_term + HB_term;
      double relaxation_time = s + 1.0;

      double eta_N = (1 / (relaxation_time + regularization)) * s;

      double eta_N_deriv;

      if(n < 2.0) {
        eta_N_deriv = ((1 / (relaxation_time + regularization)) - ((pow(1, 2.0) / (pow(relaxation_time, 2.0) + regularization)) * s)) * ((- Bi + beta * (n - 1.0) * pow(second_inv, n)) / (pow(second_inv, 2.0) + regularization)) * second_inv_deriv;
      }
      else {
        eta_N_deriv = ((1 / (relaxation_time + regularization)) - ((pow(1, 2.0) / (pow(relaxation_time, 2.0) + regularization)) * s)) * (((- Bi / (pow(second_inv, 2.0) + regularization)) + beta * (n - 1.0) * pow(second_inv, n - 2.0))) * second_inv_deriv;
      }

      // Modify log-conform parameters
      eta_p[] = s - eta_N - Wi * relaxation_time * eta_N_deriv;
      lambdav[] = Wi * relaxation_time;
    }
  }
}

// Ensure correct time for calculations
// Necessary because variable t has unreliable behavior in iterations
event update_time (i ++) {
  time_variable = t;
}

// Set mu from centered.h equal to eta_N
event update_mu (i ++) {

  double shear_rate;

  double second_inv;

  #ifdef ANALYTICAL

    #ifdef OSCILLATORY
      shear_rate = cos(time_variable / gamma_0);
    #else
      shear_rate = 1.0;
    #endif

    // Calculate second invariant for shear experiments analytically
    second_inv = sqrt(2) * fabs(shear_rate);
  
  #else

    double dux_dx = (interpolate(u.x, 0.0, 0.0) - interpolate(u.x, -1.0, 0.0));
    double dux_dy = (interpolate(u.x, 0.0, 0.0) - interpolate(u.x, 0.0, -1.0));
    double duy_dx = (interpolate(u.y, 0.0, 0.0) - interpolate(u.y, -1.0, 0.0));
    double duy_dy = (interpolate(u.y, 0.0, 0.0) - interpolate(u.y, 0.0, -1.0));

    // Calculate second invariant for 2D in general
    second_inv = sqrt(4 * pow(dux_dx, 2) + 4 * pow(duy_dy, 2) + 2 * pow(dux_dy + duy_dx, 2));
  #endif

  double yield_term = Bi / (second_inv + regularization);

  double HB_term;

  if(n < 1.0) {
    HB_term = (beta * pow(second_inv, n)) / (second_inv + regularization);
  }
  else {
    HB_term = beta * pow(second_inv, n - 1.0);
  }

  double s = yield_term + HB_term;
  double relaxation_time = s + 1.0;

  // Set mu of centered.h to correct value
  // Syntax because we use two-phase.h
  mu1 = mu2 = (1 / (relaxation_time + regularization)) * s;
}

