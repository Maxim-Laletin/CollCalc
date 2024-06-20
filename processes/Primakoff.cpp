#include <string>
#include "Model.h"

const double Mmu = 0.1056;
const double MPhoton = 0.00000001; 

// The user definitely needs to provide these four masses
double Mi = Mmu;
double Mj = MPhoton;
double Mk = Mmu;
double Mx = 0.000000001;

char s_i = 'f';
char s_j = 'b';
char s_k = 'f';

std::string pr_name = "Prim";

// Amplitude squared (without the constant part in front)
double M2(double s, double t)
{
	return t*t/(s - Mi*Mi)/(s + t - Mi*Mi);
}

