//  file: GslOde.cpp
// 
//  Definitions for the GSL versions of the Ode and Rhs C++ classes. 
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//     02/10/09  Original version using ode_test.cpp as a guide.
//
//*****************************************************************
// include files
#include <iostream>
#include <string>     // C++ strings 
#include <cstdlib>   // this has exit                                
#include <gsl/gsl_errno.h>   // header for gsl error handling
#include <gsl/gsl_matrix.h>  // header for gsl matrices
#include <gsl/gsl_odeiv.h>   // header for gsl ode solvers

#include "GslOde.h"      // include the header for these classes

//********************************************************************

// Constructor for Ode
Ode::Ode(Rhs &passed_rhs, double eps_abs, double eps_rel, std::string type) 
{
  ode_type = type;    // set the private variable for the Ode type
  ode_dim = passed_rhs.get_num_eqs();
 
  // Allocate the ode according to ode_type
  // some possibilities (see GSL manual):          
  //    gsl_odeiv_step_rk4;
  //    gsl_odeiv_step_rkf45;
  //    gsl_odeiv_step_rkck;
  //    gsl_odeiv_step_rk8pd;
  //    gsl_odeiv_step_rk4imp;
  //    gsl_odeiv_step_bsimp;  
  //    gsl_odeiv_step_gear1;
  //    gsl_odeiv_step_gear2;
  //
  if (ode_type == "rk4")
  {
    type_ptr = gsl_odeiv_step_rk4;
  }
  else if (ode_type == "rk45")
  {
    type_ptr = gsl_odeiv_step_rkf45;
  }
  else
  {
    std::cout << "Illegal ode type for GSL!" << std::endl;
    exit (1);  // time to quit!
  }

  // Allocate/initialize the stepper, the control function, and the
  //  evolution function.
  step_ptr = gsl_odeiv_step_alloc (type_ptr, ode_dim);
  control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
  evolve_ptr = gsl_odeiv_evolve_alloc (ode_dim);

  // Load values into the ode_system structure 
  ode_system.function = passed_rhs.gsl_rhs;       // functions dy[i]/dt 
  ode_system.jacobian = passed_rhs.gsl_jacobian;  //  Jacobian df[i]/dy[j] 
  ode_system.dimension = ode_dim;    // number of diffeq's 
  ode_system.params = &passed_rhs;   // pass the Rhs object
}

Ode::~Ode () // Destructor for Ode
{
  // all done; free up the gsl_odeiv stuff 
  gsl_odeiv_evolve_free (evolve_ptr);
  gsl_odeiv_control_free (control_ptr);
  gsl_odeiv_step_free (step_ptr);
}

void 
Ode::evolve (double *t, double t_next, double *h, double *y)
{
  gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                          &ode_system, t, t_next, h, y);
}

//********************************************************************

int
Rhs::gsl_rhs (double t, const double y[], double f[], void * Rhs_ptr)
{
  ((Rhs*)Rhs_ptr)->rhs (t, y, f);
  return 0;    // successful completion
}

int
Rhs::gsl_jacobian (double t, const double y[], double *dfdy,
                      double dfdt[], void * Rhs_ptr)
{
  int ode_dim = ((Rhs*)Rhs_ptr)->get_num_eqs();
  ((Rhs*)Rhs_ptr)->dfdy_mat = gsl_matrix_view_array(dfdy, ode_dim, ode_dim);
  // m_ptr points to the matrix
  ((Rhs*)Rhs_ptr)->m_ptr = &((Rhs*)Rhs_ptr)->dfdy_mat.matrix;  

  ((Rhs*)Rhs_ptr)->jacobian (t, y, dfdy, dfdt);
  return 0;   // successful completion
}

void 
Rhs::set_jacobian (int i, int j, double value)
{
  // Define the Jacobian matrix using GSL matrix routines.
  //  (see the GSL manual under "Ordinary Differential Equations") 
  gsl_matrix_set (m_ptr, i, j, value);
}

//********************************************************************
