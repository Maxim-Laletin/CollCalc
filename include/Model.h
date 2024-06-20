#ifndef __MODEL__
#define __MODEL__

#include <string>

// These values are initialized in the process file
extern double Mi, Mj, Mk, Mx; 	// masses of particles
extern char s_i, s_j, s_k; 	// types of particles

extern std::string pr_name; 	// name of the process

// This function is defined in the process file
double M2(double s, double t); 	// Amplitude squared of the process

#endif
