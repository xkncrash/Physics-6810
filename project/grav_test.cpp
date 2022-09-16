//  file: grav_test.cpp
// 
//  Simulates N-Body gravitational interactions given masses and 
//  initial values of the bodies. Numerically solves the Ode's of 
//  gravitational interaction using Gsl's Ode solver. Rhs is given 
//  by velocity and acceleration, derived from the gravitational 
//  potential energy equation. 
//
//  Programmer:  Kyle Neumann    neumann.110@osu.edu
//
//  Revision history:
//     04/16/22  original version, translated from python successfully
//     04/20/22  Implementation of GnuplotPipe
//     04/25/22  Fixed stutter due to calculation
//
//
//  Notes:
//   * Originating in a project for Theoretical Mechanics by 
//     Dr. Furnstahl of 2-body attraction where it was finished 
//     and modified to 3-body attraction for said project
//   * Modified again for a personal project to work in 
//     astronaumical units and n-bodies with a visually pleasing 
//     export gif
//   * Will be compared with the known period for Earth and 
//     analyticaly 2-body orbits to gauge error
//   * Use make -f make_grav_test to generate executable 
//
//
//  To-Do:
//   1) Comment on everything important
//   2) Fufill complete requirements of the project
//   3) Remove need for plt file or create universal plt file
//   4) Replicate HO cpp where it outputs graph without plt and 
//      in time
//   5) Generate .csv to allow use in python to utilize gif feature 
//      and better plotting software. 
//
//*****************************************************************

// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>    // C++ stringstream class (can omit iostream)
#include <cmath>
using namespace std;

#include "GslOde.h"   // Ode class for gsl (include Ode and Rhs)
#include "GnuplotPipe.h"

#define sqr(x) ((x)*(x))

// Units and Constants

const double G = 6.67408 * pow(10,-11);   // (Gravitation Const) m^3 kg^-1 s^-2
const double EM = 5.9722 * pow(10,24);    // (Earth Masses) kg 
const double SM = 1.9885 * pow(10,30);           // (Solar Masses) kg
const double yr = 365.256363004 * 86400.; // (Year) s
const double kmps = 1000.;                // (Kilometers Per Second) m/s
const double AU = 149597870691.;          // (Astronaumical Units) m
const double pi = 3.1415926535;

class Rhs_VdP : public Rhs
{
public:
  Rhs_VdP (double* M_passed, int N_passed) {N = N_passed; num_eqs = 4*N_passed; for (int i = 0; i < N_passed; i++){M[i] = M_passed[i];}};
  ~Rhs_VdP () {};
  virtual int rhs (double t, const double y[], double f[]);
  virtual int jacobian (double t, const double y[], double *dfdy, 
			double dfdt[]);
  int get_N () {return N;};
  double *get_M () {return M;}
private:
  int N;      // Number of Objects
  double M[]; // Masses
};

// function prototypes
int evolve_and_print(Ode &vdp_ode, Rhs_VdP &vdp_rhs, double* y, 
		     const double tmin, const double tmax, 
		     const double delta_t, bool fixed_CoM, bool earth_test);
double CoM_fix (double *y, double *CoM, Rhs_VdP &vdp_rhs, bool activate);


int
main ()
{
  const double eps_abs = 1.e-8;    // absolute error requested 
  const double eps_rel = 1.e-10;   // relative error requested 

  // Time parameters
  double tmin = 0.*yr;
  double tmax = 20. *yr;            
  double delta_t = 0.0001 *yr;     // Determines number of output points
  
  bool fixed_CoM = 1;              // Determines if we want CoM fixed
  bool earth_test = 0;             // Decides to run analytical comparison to Sun-Earth System
  bool three_body = 0;             // Preset up variables for 3-body orbit

  // Initial and constant variables
  // Interesting 4-Body system
  
  // Masses
  double mass[] = {1.* SM,1.* SM,0.5 * SM, 1000. * EM};
      
  // Initial radius from an irrelevant origin
  double r_0[] = {1.*AU,5.*AU, 3.*AU, 6. *AU}; 

  // Initial phase in relation to origin             
  double theta_0[] = {0, pi, .5*pi, 1.9*pi};   
 
  // Initial tangential              
  double v_0[] = {6.*kmps,8.*kmps, 20.*kmps, 12.*kmps};  
  
  int N = sizeof(mass)/sizeof(mass[0]);            // Number of objects

  // Verify all of our arrays contain the same number of objects
  if (sizeof(r_0)/sizeof(r_0[0]) != N)
    {
      cout << "Initial variable length doesn't agree with N" << endl;
      return 0;
    }
  if (sizeof(v_0)/sizeof(v_0[0]) != N)
    {
      cout << "Initial variable length doesn't agree with N" << endl;
      return 0;
    }
  if (sizeof(theta_0)/sizeof(theta_0[0]) != N)
    {
      cout << "Initial variable length doesn't agree with N" << endl;
      return 0;
    }

  // If earth_test is active, this will overwrite the arrays into only 
  // the Sun-Earth system (approximately) to test the accuracy of the 
  // minute 2-body system.
  if (earth_test)
    {
      // Set some given parameters of the system
      mass[0] = 1. * SM;
      mass[1] = 1. * EM;
      r_0[0] = mass[1]*AU/mass[0];     // Uses center of mass rules
      r_0[1] = AU;

      double tperiod = sqrt(4.*sqr(pi) * AU * sqr(AU)/(G*mass[0]));   // Kepler's 3rd Law

      v_0[0] = 2.*pi*r_0[0]/tperiod;   // Distance/time
      v_0[1] = 2.*pi*r_0[1]/tperiod;
      theta_0[0] = 0.;
      theta_0[1] = pi;
      N = 2;             // Crucial to prevent any other variables of other objects from being read

      tmax = 1.*tperiod;     // Test how the system reacts over time
    }
  else if (three_body)
    {
      // Sets up initial parameters for a 3-body system 
      mass[0] = 2. * SM;
      mass[1] = .4 * SM;
      mass[2] = .2 * SM;
      r_0[0] = 2*mass[1]*AU/mass[0];     // Uses center of mass rules
      r_0[1] = 2*AU;
      r_0[2] = 6.* AU;
      v_0[0] = 10. *kmps;   // Distance/time
      v_0[1] = 20. * kmps;
      v_0[2] = 17. * kmps;
      theta_0[0] = 0.;
      theta_0[1] = pi;
      theta_0[2] = 1.4*pi;
      N = 3;   
    }

  // Converts polar coordinates into cartesian
  double y[4*N];       // Array of all positional and velocity variables
  for (int i = 0; i < N; i++)
    {
      y[4*i] = r_0[i]*cos(theta_0[i]);
      y[4*i+1] = -v_0[i]*sin(theta_0[i]);
      y[4*i+2] = r_0[i]*sin(theta_0[i]);
      y[4*i+3]= v_0[i]*cos(theta_0[i]);
    }

  // Creates the Rhs_VdP object
  Rhs_VdP vdp_rhs (mass, N);

  // Creates the Ode object from GslOde.cpp
  Ode vdp_ode (vdp_rhs, eps_abs, eps_rel, "rk45");

  // Uses the function to evolve our variables in time
  evolve_and_print(vdp_ode, vdp_rhs, y, 
		   tmin, tmax, delta_t, fixed_CoM, earth_test);

  return 0;
}

//*************************************************************

//********************evolve_and_print*************************
//
// Evolves array y in time using Gsl's Ode solver. 
// Outputs the resulting new coordinates into y 
// and output file.
//
//*************************************************************
int
evolve_and_print(Ode &vdp_ode, Rhs_VdP &vdp_rhs, double* y, 
		     const double tmin, const double tmax, 
		 const double delta_t, bool fixed_CoM, bool earth_test)
{
  // Grabs variables from vdp_rhs object
  int N = vdp_rhs.get_N();
  double *M = vdp_rhs.get_M();
  double CoM[2] = {0};

  bool run_gnuplot = 1;    // Determines if the gnuplot will run, on unless something breaks

  // declare a GnuplotPipe object and set some properties
  int pts = 100+50*(tmax-tmin)/yr;
  double delta_t_plot = 1.*(tmax-tmin)/pts;
  double t_next_plot = 0.;

  int total_steps = ((tmax-tmin))/delta_t_plot;
  double pos[N][2][total_steps+1];
  int t_step = 0;

  ostringstream my_gnutitle;
  my_gnutitle << N << "-Body Orbit Over " << (tmax-tmin)/yr << " Years"; 
  GnuplotPipe myPipe (N);
  myPipe.set_title (my_gnutitle.str());
  myPipe.set_xlabel ("x");
  myPipe.set_ylabel ("y");

  // Used to determine a good plot size to keep scaling
  double max_R = 0, r = 0;
  for (int i = 0; i < N; i++)
    {
      r = sqrt(sqr(y[4*i])+sqr(y[4*i+2]));
      
      if (r/AU > max_R){max_R = r/AU;}
    }
  myPipe.set_xmin (-1.5 * max_R);
  myPipe.set_xmax (1.5 * max_R);
  myPipe.set_ymin (-1.5 * max_R);
  myPipe.set_ymax (1.5 * max_R);

  myPipe.set_delay(5000000/pts);
  /*
  if ((100 * delta_t/yr) < 1)
    {
  myPipe.set_delay(round(100 * delta_t/yr));
    }
  else {myPipe.set_delay(0);}
  */
  if (run_gnuplot)
    {
      myPipe.init (M);  // Start up piping to gnuplot
    }
  
  // Sets initial time 
  double t = tmin;

  // Generates output file given number of objects
  ostringstream my_stringstream;
  my_stringstream << N << "_orbital.dat";
  ofstream my_out;
  my_out.open(my_stringstream.str().c_str());

  my_out << "# Running gravitational simulaition with N = " << N 
	 << " bodies" << endl;
  my_out << "#" ;
  for (int i = 0; i < N; i++)
    {
      my_out << "  m_" << i+1 << " = " << M[i];   
    }
  my_out << endl;
  my_out << "# t";
  for (int i = 0; i < N; i++)
    {
      my_out << "  x_" << i+1 << "  y_" << i+1; 
    }

  my_out << endl;
  my_out << scientific << setprecision (5) << setw (12) << t/yr << " ";
  
  // Uses CoM function to prevent the motion of unstable systems in 
  // terms of plot if fixed_CoM is true
  
  CoM_fix(y, CoM, vdp_rhs, fixed_CoM);
    

  // Initial coordinates with CoM centering
  for (int i = 0; i < N; i++)
    {
      my_out << setw (12) << (y[4*i]-CoM[0])/AU << " " 
	     << setw (12) << (y[4*i+2]-CoM[1])/AU << " ";
      pos[i][0][t_step] = (y[4*i]-CoM[0])/AU;
      pos[i][1][t_step] = (y[4*i+2]-CoM[1])/AU;
      /*
      if (run_gnuplot)
	{
	  myPipe.plot ((y[4*i]-CoM[0])/AU,(y[4*i+2]-CoM[1])/AU,i);
	  if (i == N-1) {t_next_plot += delta_t_plot;}
	}
      */
    }
  t_step++;
  my_out << endl;

  double h = 1e-6;            // Step size Gsl takes to evolve t, not output

  // Separates set up if we want a fixed CoM. 
  // More efficient than looped if statements
   
  for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
      while (t < t_next)  // evolve from t to t_next 
	{
	  vdp_ode.evolve ( &t, t_next, &h, y );    
	}
      CoM_fix ( y, CoM, vdp_rhs, fixed_CoM);

      my_out << scientific << setprecision (5) << setw (12) << t /yr<< " ";

      for (int i = 0; i < N; i++)
	{
	  my_out << setw (12) << (y[4*i]-CoM[0])/AU << " " 
		 << setw (12) << (y[4*i+2]-CoM[1])/AU << " ";

	  // Prevents too many points from slowing down system
	  if (run_gnuplot && t >= t_next_plot)
	    {
	      pos[i][0][t_step] = (y[4*i]-CoM[0])/AU;
	      pos[i][1][t_step] = (y[4*i+2]-CoM[1])/AU;
	      if (i == N-1){t_next_plot += delta_t_plot; t_step++;}
	    }
	  
	}

      my_out << endl;
    }
  if (run_gnuplot)
    {
      for (int j = 0; j < total_steps+1; j++)
	{
	  for (int i = 0; i < N; i++)
	    {
	      myPipe.plot (pos[i][0][j],pos[i][1][j],i);
	    }
	}
    }
  /*
  else  
    {
      for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
	{
	  while (t < t_next)  // evolve from t to t_next 
	    {
	      vdp_ode.evolve ( &t, t_next, &h, y );
	    }

	  my_out << scientific << setprecision (5) << setw (12) 
		 << t /yr<< " ";

	  for (int i = 0; i < N; i++)
	    {
	      my_out << setw (12) << y[4*i]/AU << " " << setw (12) 
		     << y[4*i+2]/AU << " ";
	      
	      // Prevents too many points from clogging system
	      if (run_gnuplot && t >= t_next_plot)
		{
		  myPipe.plot ((y[4*i])/AU,(y[4*i+1])/AU,i);
		  if (i == N-1) {t_next_plot += delta_t_plot;}
		}
	    }
	  my_out << endl;
	}
    }
  */
  // Tests the final point of the Earth to see how close it came to 
  // making 1 full orbit in the analytical time calculated
  if (earth_test)
    {
      double coord[2] = {-AU,0.};    // Analytical final position

      double earth[2] = {y[4],y[6]};    // Final numerical position

      double relerror = sqrt(sqr(coord[0]-earth[0])+sqr(coord[1]-earth[1])
			     )/(-coord[0]);    // Relative error of position
      cout << endl << setprecision (5) << fixed;
      cout << "Comparing to an Sun-Earth system, our relative error is " 
	   << relerror * 100. << "%" << endl;
      cout << "Final position of (" << y[4]/AU << " AU, " << y[6]/AU 
	   << " AU) " << "at time " << t/yr << " years" << endl;   
      cout << "Final step size used to find numerical positions was " 
	   << h << endl << endl;
    }
    
  if (run_gnuplot){myPipe.finish();}
  my_out.close();

  cout << "All files outputted" << endl;

  return (0);
}

//*****************************************************************

//**************************CoM_fix*********************************
//
// Determines the center of mass (CoM) of the systm, so it 
// can be subtracted from the output position, fixing the 
// CoM at the origin.
//
//*******************************************************************
double 
CoM_fix (double *y, double * CoM, Rhs_VdP &vdp_rhs, bool activate)
{
  // Grabs variables from object
  if(activate)
    {
      int N = vdp_rhs.get_N();
      double *M = vdp_rhs.get_M(); 
      double tot_M = 0., x_c = 0., y_c = 0.;         // Total mass, x CoM, y CoM
      for (int n = 0; n < N; n++)
	{
	  tot_M += M[n];
	  x_c += M[n] * y[4*n];
	  y_c += M[n] * y[4*n + 2];
	}
      x_c /= tot_M;
      y_c /= tot_M;

      CoM[0] = x_c;
      CoM[1] = y_c;
    }

  return (0);
}
//************************************************************

//************************************************************

//*************************rhs*******************************
//
// Right hand side equation of the ODE to be 
// numerically solved.
//
//************************************************************
int
Rhs_VdP::rhs (double, const double y[], double f[])
{
  //double f [4*N] = {0};

  // loops for each object and every object it interacts with
  for (int a = 0; a < N; a++)
    {
      double P = 0., Q = 0.;

      for (int b = 0; b < N; b++)
	{
	  if (a!=b)  // Ignores when a = b as variables would be zero or worse
	    {
	      //Distance separation between object a and b cubed
	      double d_sqr = pow(sqr(y[4*a]-y[4*b])+sqr(y[4*a+2]-y[4*b+2])
				 ,3./2.); 
	      P += M[b]*(y[4*a]-y[4*b])/d_sqr;      // x acceleration term
	      Q += M[b]*(y[4*a+2]-y[4*b+2])/d_sqr;  // y acceleration term
	    }
	}
      // Rhs array
      f[4*a] = (y[4*a+1]);
      f[4*a+1] = -G*P;
      f[4*a+2] = (y[4*a+3]);
      f[4*a+3] = -G*Q;
    }

  return 0;
}

int 
Rhs_VdP::jacobian (double , const double y[], double *,	double dfdt[])
{
  // Ignore completely
  set_jacobian (0,0,y[0]);
  for (int i = 0; i < 4*N; i++)
    {
dfdt[i] = 0.0;
    }
  return (0);
}
