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

#define sqr(x) ((x)*(x))

// Units and Constants

const double G = 6.67408 * pow(10,-11);   // (Gravitation Const) m^3 kg^-1 s^-2
const double EM = 5.9722 * pow(10,24);    // (Earth Masses) kg 
const double SM = 332900. * EM;           // (Solar Masses) kg
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
		     const double delta_t, bool fixed_CoM);
double CoM_fix (double *y, double *CoM, Rhs_VdP &vdp_rhs);

int
main ()
{
  const double eps_abs = 1.e-8;    // absolute error requested 
  const double eps_rel = 1.e-10;   // relative error requested 

  // Time parameters
  double tmin = 0.*yr;
  double tmax = 10. *yr;
  double delta_t = 0.0001 *yr;     // Determines how many output times

  bool fixed_CoM = 1;              // Determines if we want CoM fixed

  // Initial and constant variables
  double mass[] = {1.* SM,1.* SM,0.5 * SM, 1000. * EM};      // Masses
  double r_0[] = {1.*AU,0.5*AU, 3.*AU, 5. *AU};              // Initial radius from irrelevant origin
  double theta_0[] = {0, pi, .5*pi, .9*pi};                  // Initial phase in relation to origin
  double v_0[] = {14.9*kmps,14.9*kmps, 30.*kmps, 20.*kmps};  // Initial tangential

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
		   tmin, tmax, delta_t, fixed_CoM);

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
		 const double delta_t, bool fixed_CoM)
{
  // Grabs variables from vdp_rhs object
  int N = vdp_rhs.get_N();
  double *M = vdp_rhs.get_M();
  double CoM[2] = {0};

  /*  double y[4*N];    // current solution vector
  for (int i = 0; i < N; i++)
    {
      y[4*i] = x0[i];
      y[4*i+1] = x_dot_0[i];
      y[4*i+2] = y0[i];
      y[4*i+3] = y_dot_0[i];
      }*/

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
  if (fixed_CoM)
    { 
      CoM_fix(y, CoM, vdp_rhs);
    }  

  // Initial coordinates with CoM centering
  for (int i = 0; i < N; i++)
    {
      my_out << setw (12) << (y[4*i]-CoM[0])/AU << " " << setw (12) << (y[4*i+2]-CoM[1])/AU << " ";
    }

  my_out << endl;

  double h = 1e-6;            // Step size Gsl takes to evolve t, not output

  // Separates set up if we want a fixe CoM. 
  // More efficient than looped if statements
  if (fixed_CoM)  
    { 
      for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
	{
	  while (t < t_next)  // evolve from t to t_next 
	    {
	      vdp_ode.evolve ( &t, t_next, &h, y );    
	    }
	  CoM_fix ( y, CoM, vdp_rhs);

	  my_out << scientific << setprecision (5) << setw (12) << t /yr<< " ";

	  for (int i = 0; i < N; i++)
	    {
	      my_out << setw (12) << (y[4*i]-CoM[0])/AU << " " << setw (12) << (y[4*i+2]-CoM[1])/AU << " ";
	    }

	  my_out << endl;
	}
    }
  else  
    {
      for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
	{
	  while (t < t_next)  // evolve from t to t_next 
	    {
	      vdp_ode.evolve ( &t, t_next, &h, y );
	    }

	  my_out << scientific << setprecision (5) << setw (12) << t /yr<< " ";

	  for (int i = 0; i < N; i++)
	    {
	      my_out << setw (12) << y[4*i]/AU << " " << setw (12) << y[4*i+2]/AU << " ";
	    }

	  my_out << endl;
	}
    }

  my_out.close();

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
CoM_fix (double *y, double * CoM, Rhs_VdP &vdp_rhs)
{
  // Grabs variables from object
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
