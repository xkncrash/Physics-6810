//  file: GslOde.h
//
//  Header file for the GSL versions of the Ode and Rhs C++ classes. 
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02/10/09  original version
//
//  Notes:
//   * We encapsulate GSL ode functions in this class, based on the
//      documentation for the GSL library under "Differential Equations".
//   * The GSL header file is included in GslOde.cpp and
//      anywhere we want to create Ode and Rhs objects.
//   * See the GSL documentation for error handling details.
//
//*****************************************************************
// The ifndef/define macro ensures that the header is only included once
#ifndef GSLODE_H
#define GSLODE_H

// include files
#include <string>     // C++ strings 
                                
#include <gsl/gsl_odeiv.h>   // header for gsl ode solvers
#include <gsl/gsl_matrix.h>  // header for gsl matrices

class Rhs
{
  public:
    Rhs () {};
    virtual ~Rhs () {};
    virtual int rhs (double, const double *, double *) {return (1);};
    virtual int jacobian (double, const double *, double *, 
                          double *) {return (1);};
    static int gsl_rhs (double t, const double y[], double f[], void *);
    static int gsl_jacobian (double t, const double y[], double *dfdy,
                             double dfdt[], void *);
    int get_num_eqs () {return num_eqs;};
  
  protected:  
    int num_eqs;                  
    void set_jacobian (int i, int j, double value);
    gsl_matrix_view dfdy_mat;
    gsl_matrix *m_ptr;  
                        
  private:
};

class Ode
{ 
  public:
    Ode (Rhs &this_rhs, double eps_abs, double eps_rel, std::string type);  
    ~Ode ();  
    void evolve (double *t, double t_next, double *h, double y[]);

  private: 
    int ode_dim;                // number of equations
    std::string ode_type;       // type of Ode 
    const gsl_odeiv_step_type *type_ptr;
    // The GSL steeper, control function, and the evolution function.
    gsl_odeiv_step *step_ptr;
    gsl_odeiv_control *control_ptr;
    gsl_odeiv_evolve *evolve_ptr;

    gsl_odeiv_system ode_system;  // structure with the rhs function, etc. 
};

#endif
