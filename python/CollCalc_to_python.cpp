#include <CollCalc_to_python.h>

// ************************  Main file  ************************ 

// Assign +/- or 0 depending on the statistics of the particle
double assign_sign(char type)
{
	switch(type){
		case 'f': // Fermi-Dirac
			return 1.0;
			break;
 		case 'b': // Bose-Einstein
			return -1.0;
			break;
		case 'm': // Maxwell-Boltzmann
			return 0.0;
			break;
		default :
			return 0.0;
	}
}

	// Find the mass of the heaviest particle (to define x)
	double Mmax = std::max(std::max(Mi,Mj),std::max(Mk,Mx));

	// Assign the signs to the distribution functions (DF)
	double sign_fi = assign_sign(s_i);
	double sign_fj = assign_sign(s_j);
	double sign_fk = assign_sign(s_k);

	// Desired relative accuracy of the integration
	double rel_acc_ek, rel_acc_cos_s;
	double rel_acc_cos_t = 0.1;
	double rel_acc_cos_phi = 0.1;

	// Size of memory allocated for each integration
	int mem_alloc_ek = 1000;
	int mem_alloc_cos_s = 1000;
	int mem_alloc_cos_t = 1000;
	int mem_alloc_cos_phi = 1000;

	// Keys that correspond to different number of nodes in the GK quadrature
	int gsl_GK_key_ek = 2;
	int gsl_GK_key_cos_s = 1;
	int gsl_GK_key_cos_t = 1;
	int gsl_GK_key_cos_phi = 1;


// ========== Functions to be imported into Python ==========

// Collision integral function for one particle with unknown distribution
double CoAnnihilation(double x, double q, double rel_acc)
{
	
	double rel_err, result; // relative error of the integration and result

	// Desired relative accuracy of the integration
	rel_acc_ek = rel_acc;
	rel_acc_cos_s = 0.1;

	// The function-pointer points to the function M2(s,t) defined in the model.cpp file
	Msquared = M2;

	result = integral_ek(x,q,rel_err);

	return result;
}


// Collision integral function for two particles with the same unknown distribution
double Annihilation(double x, double q, double p, double rel_acc)
{

	double rel_err, result; // relative error of the integration and result

	// Desired relative accuracy of the integration
	rel_acc_cos_s = rel_acc;

	// The function-pointer points to the function M2(s,t) defined in the model.cpp file
	Msquared = M2;

	result = integral_cos_s(x,q,p,rel_err);

	return result;
}
