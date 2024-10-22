// Example process file for the annihilation of muons into photon and axion
#include "Model.h"

// **********************  Process file  **********************

// Muon mass
const double Mmu = 0.1056; // [GeV]
// Photon mass
const double MPhoton = 0.00001; // [GeV] non-zero mass to avoid potential errors

// Masses of the particles
double Mi = Mmu;
double Mj = Mmu;
double Mk = MPhoton;
double Mx = 0.00001; // axion mass [GeV]

// Types of particles
char s_i = 'f';
char s_j = 'f';
char s_k = 'b';

// Name of the process
std::string pr_name = "Ann";

// Amplitude squared (without the constant part in front)
double M2(double s, double t)
{
	return s*s/(Mmu*Mmu - t)/(s + t - Mmu*Mmu); // [1]
}
