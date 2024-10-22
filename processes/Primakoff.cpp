// Example process file for the Primakoff scattering of axion on muon
#include <string>
#include "Model.h"

const double Mmu = 0.1056;
const double MPhoton = 0.00000001; // non-zero mass to avoid potential errors

double Mi = Mmu;
double Mj = MPhoton;
double Mk = Mmu;
double Mx = 0.000000001;

char s_i = 'f';
char s_j = 'b';
char s_k = 'f';

// Process name used to create the file with the results
std::string pr_name = "Prim";

// Amplitude squared (without the constant part in front)
double M2(double s, double t)
{
	return t*t/(s - Mi*Mi)/(s + t - Mi*Mi);
}

