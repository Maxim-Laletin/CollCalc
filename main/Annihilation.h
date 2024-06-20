#ifndef __ANNIHILATION__
#define __ANNIHILATION__

#include <vector> 
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "Functions.h"
#include "Model.h"

// These values are initialized in the main file
extern double Mmax, sign_fi, sign_fj, sign_fk;
extern double rel_acc_ek, rel_acc_cos_s, rel_acc_cos_t, rel_acc_cos_phi;
extern int mem_alloc_ek, mem_alloc_cos_s, mem_alloc_cos_t, mem_alloc_cos_phi; 
extern int gsl_GK_key_ek, gsl_GK_key_cos_s, gsl_GK_key_cos_t, gsl_GK_key_cos_phi;

// This function pointer is initialized in the main file
inline double (*Msquared)(double,double);
// (the function that it points to is defined in the process file)


// The base integral for annihilation (over cos_s)
double integral_cos_s (double x, double q, double p, double &relerror);

// The base integral for co-annihilation (over ek)
double integral_ek (double x, double q, double &relerror);

// Integrand functions (same structure to construct GSL functions)
double integrand_ek (double ek, void * params);

double integrand_cos_s (double cos_s, void * params);

double integrand_cos_t (double cos_t, void * params);

double integrand_cosphi (double cosphi, void * params);

#endif
