#include "Model.h"

// **********************  Process file  **********************

// Constant values of masses (MAYBE TRANSFERRED TO A DIFFERENT FILE?)
const double Mmu = 0.1056; // [GeV]
const double MPhoton = 0.00001; // [GeV]

// Masses of the particles
double Mi = Mmu;
double Mj = Mmu;
double Mk = MPhoton;
double Mx = 0.00001; // DM particle mass [GeV]

// Types of particles
char s_i = 'f';
char s_j = 'f';
char s_k = 'b';

// Name of the process
std::string pr_name = "Ann";

// Amplitude squared 	(without the constant part in front)
double M2(double s, double t)
{
	return s*s/(Mmu*Mmu - t)/(s + t - Mmu*Mmu); // [1]
}
