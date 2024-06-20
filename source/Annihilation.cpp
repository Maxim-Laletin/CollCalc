#include "Annihilation.h"

// **********************  Annihilations and co-annihilations **********************

//   Scheme of the process and definitions
/*
  (x)  q 			  			p_j   (j)
				\						 /				
				 \		   		/
				  *>----->*
				 /		   		\
				/						 \
  (k)  p 			  			p_i   (i)
*/
// (it's just a sketchy representation of the process kinematics, not a Feynman diagram of the s-channel reaction)


// Variable to dump the error from all gsl integrations
static double error_bulk; 


// **********************  Annihilations  **********************

// For annihilation m_k = m_x, so both k and x represent a particle with an unknown PDF
// and i and j are some particles with an equilibrium PDF


// Keep the integration parameters initialized in the root integral
struct integration_parameters
{
	gsl_integration_workspace * w_cos_s;
	gsl_integration_workspace * w_cos_t;
	gsl_integration_workspace * w_cos_phi;
	int gsl_GK_key_cos_s; 
	int gsl_GK_key_cos_t; 
	int gsl_GK_key_cos_phi;
};

// Structure of parameters for the integrand functions
struct general_parameters
{
	std::vector<double> pars;
	struct integration_parameters * ipars; // pointer to a structure of integration parameters
};


// The base integral for annihilation (over cos_s)
double integral_cos_s (double x, double q, double p, double &relerror)
{

	// Inversed temperature
	double T = Mmax/x;
	// Reduced masses of particles (m/T)
	double mi = Mi/T;
	double mj = Mj/T;
	double mk = Mk/T;
	double mx = Mx/T;

	// Energy of x particle
	double ex = sqrt(mx*mx + q*q);
	// Energy of k particle
	double ek = sqrt(mk*mk + p*p);

	// The sum of the energies
	double esum = ek+ex;

	// The part of s variable that is independent of cos_s (see the definition of s in the integrand below)
	double sbase = mx*mx + mk*mk + 2*ex*ek;

	double cos_s_max = fmin( 1.0, (sbase - pow(mi+mj,2))/(2*q*p) );
	// cos_s_max is obtained from the condition that s(cos_s) > (mi+mj)^2

	// CHECK ALL THE PROPER INTEGRATION LIMITS HERE

	double cos_s_min = -1.0; // CHECK AND CHANGE !!!
	

	double integral_cos_s_res = 0.0; // store the result of integration
	double quad_error = 0.0; // store the absolute error of the integration
	relerror = 0.0; // store the relative error of the integration

	if (cos_s_min < cos_s_max) // CHECK AND CHANGE !!!
	{


		int gsl_GK_key_cos_s = 1;
		int gsl_GK_key_cos_t = 1;
		int gsl_GK_key_cos_phi = 1;
	
		// // (WE NEED TO IMPOSE THAT THE USER CAN SWITCH THIS PROCEDURE OFF AND USE THE KEYS THAT THEY PROVIDED)
		// // Regulate the keys to enhance the speed of the integral calculation
		// if (ek_max/ek_min >= 1E+4)
		// {
		// 	gsl_GK_key_ek = 5;
		// 	gsl_GK_key_cos_s = 3;
		// } 
		// else if (ek_max/ek_min >= 1000) 
		// {
		// 	gsl_GK_key_ek = 4;
		// 	gsl_GK_key_cos_s = 3;
		// }
		// else if (ek_max/ek_min > 100)
		// {
		// 	gsl_GK_key_ek = 3;
		// 	gsl_GK_key_cos_s = 2;
		// }

		
		// Allocate memory for integration (mem_alloc is the size)
		gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (mem_alloc_cos_s);
		gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (mem_alloc_cos_t);
		gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (mem_alloc_cos_phi);
		// The GSL manual advises to allocate a separate workspace for each integration


		// Vector of physical parameters
		std::vector<double> cos_s_pars = {mi,mj,mk,mx,T*T,q,ex,p,ek,sbase,esum};

		// Parameters for numerical integration
		struct integration_parameters ipars = {w1, w2, w3, gsl_GK_key_cos_s, gsl_GK_key_cos_t, gsl_GK_key_cos_phi};

		// Creating the array of general parameters for nested integrations (pars + ipars)
		struct general_parameters g_p = {cos_s_pars, &ipars};

		gsl_function F_cos_s;
		F_cos_s.function = &integrand_cos_s;
		F_cos_s.params = &g_p; // pass the renewed structure of parameters

		gsl_integration_qag (&F_cos_s, cos_s_min, cos_s_max, 0, rel_acc_cos_s, mem_alloc_cos_s, gsl_GK_key_cos_s, w1, &integral_cos_s_res, &error_bulk); 


		// Free the memory allocated for integration
		gsl_integration_workspace_free(w1);
		gsl_integration_workspace_free(w2);
		gsl_integration_workspace_free(w3);

		relerror = quad_error/integral_cos_s_res; // relative error of the integral calculation
		// [PROBLEM CAN APPEAR IF INTEGRAL = 0]
	}

	return T*p*integral_cos_s_res/(4*pow(2*M_PI,4)); 
	
}


// **********************  Co-annihilations  **********************

// For co-annihilation only x represents a particle with an unknown PDF,
// while i, j adn k are some particles with an equilibrium PDF


// Find the (allowed) top limit of Ek
double find_ek_max(double x, double q, double ek_min, double mi, double mj, double mk, double mx)
{
	double ex = sqrt(mx*mx+q*q);
	double ek = ek_min;

	// in case ei_plasma is close to 0 and exp -> Inf
	const double machine_zero = 1E-10;

	// The middle of the cos_s range (or the maximal value, if it is negative)
	double cos_s_mid = fmin(machine_zero, 0.5*(mx*mx + mk*mk + 2*ex*ek - pow(mi+mj,2))/q/sqrt(ek*ek-mk*mk) + machine_zero );
	double s = mx*mx + mk*mk + 2*ex*ek - 2*q*sqrt(ek*ek-mk*mk)*cos_s_mid;

	double ei_plasma = e_CM(s,mi,mj)*( (ek+ex)/sqrt(s) );

	// the product of the PDFs at minimal ek
	double pdf_prod_max = 1.0/( (exp(ei_plasma)+sign_fi)*(exp(ei_plasma)+sign_fj) );
	// 
	const double max_exp_index = 15.0;
	double ek_max = max_exp_index - log(pdf_prod_max) - ex; 

	if (ek_max > 700.0 - ex) 
	{
		ek_max = 700.0 - ex;
	}

	return ek_max;
}



// The base integral for co-annihilation (over Ek)
double integral_ek (double x, double q, double &relerror)
{
	// Inversed temperature
	double T = Mmax/x;
	// Reduced masses of particles (m/T)
	double mi = Mi/T;
	double mj = Mj/T;
	double mk = Mk/T;
	double mx = Mx/T;

	// Energy of x particle
	double ex = sqrt(mx*mx + q*q);

	// Minimal energy of k particle (depends on the energy of X particle)
	// Derived from the condition that s_max > (mi+mj)^2
	double ek_min = 0.5*( ex*(pow(mi+mj,2) - mk*mk - mx*mx) - q*sqrt((pow(mx+mi+mj,2) - mk*mk)*(pow(mx-mi-mj,2) - mk*mk)) )/(mx*mx); 

	if (ek_min > 0.5*(pow(mi+mj,2) - mx*mx - mk*mk)/ex)
	{
		ek_min = mk; // if s_max is always > (mi+mj)^2
	}

	double integral_ek_res = 0.0; // store the result of integration
	double quad_error = 0.0; // store the absolute error of the integration
	relerror = 0.0; // store the relative error of the integration

	if (ek_min < 700.0 - ex)
	{

		// Generally, should be Inf
		//double ek_max = 20*ek_min;

		double ek_max = find_ek_max(x,q,ek_min,mi,mj,mk,mx);

		int gsl_GK_key_ek = 2;
		int gsl_GK_key_cos_s = 1;
		int gsl_GK_key_cos_t = 1;
		int gsl_GK_key_cos_phi = 1;
	
		// (WE NEED TO IMPOSE THAT THE USER CAN SWITCH THIS PROCEDURE OFF AND USE THE KEYS THAT THEY PROVIDED)
		// Regulate the keys to enhance the speed of the integral calculation
		if (ek_max/ek_min >= 1E+4)
		{
			gsl_GK_key_ek = 5;
			gsl_GK_key_cos_s = 3;
		} 
		else if (ek_max/ek_min >= 1000) 
		{
			gsl_GK_key_ek = 4;
			gsl_GK_key_cos_s = 3;
		}
		else if (ek_max/ek_min > 100)
		{
			gsl_GK_key_ek = 3;
			gsl_GK_key_cos_s = 2;
		}

		
		// Allocate memory for integration (mem_alloc is the size)
		gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (mem_alloc_ek);
		gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (mem_alloc_cos_s);
		gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (mem_alloc_cos_t);
		gsl_integration_workspace * w4 = gsl_integration_workspace_alloc (mem_alloc_cos_phi);
		// The GSL manual advises to allocate a separate workspace for each integration

	
		// Vector of physical parameters
		std::vector<double> ek_pars = {mi,mj,mk,mx,T*T,q,ex};

		// Parameters for numerical integration
		struct integration_parameters ipars = {w2, w3, w4, gsl_GK_key_cos_s, gsl_GK_key_cos_t, gsl_GK_key_cos_phi};

		// Creating the array of general parameters for nested integrations (pars + ipars)
		struct general_parameters g_p = {ek_pars, &ipars};

		// Initialize the ek integrand functions
		gsl_function F_ek;
		F_ek.function = &integrand_ek;
		
		// Passing the parameters to the integrand function
		F_ek.params = &g_p;

		// Integrating using the GSL adaptive GK quadrature 
		gsl_integration_qag (&F_ek, ek_min, ek_max, 0, rel_acc_ek, mem_alloc_ek, gsl_GK_key_ek, w1, &integral_ek_res, &quad_error);

		// Free the memory allocated for integration
		gsl_integration_workspace_free(w1);
		gsl_integration_workspace_free(w2);
		gsl_integration_workspace_free(w3);
		gsl_integration_workspace_free(w4);

		relerror = quad_error/integral_ek_res; // relative error of the integral calculation
		// [PROBLEM CAN APPEAR IF INTEGRAL = 0]
	}

	return (T*T)*integral_ek_res/(4*pow(2*M_PI,4)); // T^2 times the integral
}

// **********************  Integrals for both annihilations and co-annihilations  **********************


//   (the arrow shows the direction in which the integrals are nested - 
//    the integral on the bottom is executed the most number of times)
/*
									  _
								   | |
								  _| |_
								  \   /
								   \ /
								    '
*/

// ek is the energy of k particle
double integrand_ek (double ek, void * params)
{

	// Dereference the pointer to the integrand parameters
	struct general_parameters g_p = * (struct general_parameters *) params;
	//std::vector<double> pars = * (struct general_parameters *) params.pars;
	g_p.pars.reserve(11); // reserve more space for a vector of double-valued parameters

	// pars = {mi,mj,mk,mx,T^2,q,ex};

	// Create references to numerical parameters to maintain readability
	const double &mi = g_p.pars[0];
	const double &mj = g_p.pars[1];
	const double &mk = g_p.pars[2];
	const double &mx = g_p.pars[3];
	const double &q = g_p.pars[5];
	const double &ex = g_p.pars[6];


	// Momentum of p particle
	double p = sqrt(ek*ek - mk*mk);

	// The sum of the energies
	double esum = ek+ex;

	// The part of s variable that is independent of cos_s (see the definition of s in the integrand below)
	double sbase = mx*mx + mk*mk + 2*ex*ek;

	double cos_s_max = fmin( 1.0, (sbase - pow(mi+mj,2))/(2*q*p) );
	// cos_s_max is obtained from the condition that s(cos_s) > (mi+mj)^2

	// Add new values to the vector of physical parameters
	g_p.pars.insert(g_p.pars.end(), {p,ek,sbase,esum});

	double integral_cos_s;

	gsl_function F_cos_s;
	F_cos_s.function = &integrand_cos_s;
	F_cos_s.params = &g_p; // pass the renewed structure of parameters

	gsl_integration_qag (&F_cos_s, -1.0, cos_s_max, 0, rel_acc_cos_s, mem_alloc_cos_s, g_p.ipars->gsl_GK_key_cos_s, g_p.ipars->w_cos_s, &integral_cos_s, &error_bulk); 


	return p*(1.0 - sign_fk/(exp(ek)+sign_fk) )*integral_cos_s;
}

/*
									  _
								   | |
								  _| |_
								  \   /
								   \ /
								    '
*/


// cos_s is the (cosine of the) angle between p_x and p_k (s invariant angle)
double integrand_cos_s (double cos_s, void * params)
{

	struct general_parameters g_p = * (struct general_parameters *) params;
	g_p.pars.reserve(19);

	// pars = {mi,mj,mk,mx,T^2,q,ex,p,ek,sbase,esum};

	const double &mi = g_p.pars[0];
	const double &mj = g_p.pars[1];
	const double &mk = g_p.pars[2];
	const double &mx = g_p.pars[3];
	const double &q = g_p.pars[5];
	const double &ex = g_p.pars[6];
	const double &p = g_p.pars[7];
	const double &ek = g_p.pars[8];
	const double &sbase = g_p.pars[9];
	const double &esum = g_p.pars[10];


	// Mandelstam s variable (in plasma frame)
	double s = sbase - 2*q*p*cos_s;

	// Momenta in the CM frame
	double q_CM = p_CM(s,mk,mx); 
	double pi_CM = p_CM(s,mi,mj);

	// velocity of the CM frame
	//double v_CM = sqrt(1.0 - s/pow(ek+ex,2));
	// Lorenz factor
	//double y = 1/sqrt(1 - v_CM*v_CM);
	double y = esum/sqrt(s);

	// The scalar product of q and v_CM
	double q_vCM = (q*q + q*p*cos_s)/esum;

	// the angle between p_x and p_x_CM
	double cos_xCM = (q*q + y*q_vCM*(y*q_vCM/(y+1) - ex))/(q*q_CM);

	// The part of t variable that is independent of cos_t (see the definition of t in the integrand below)
	double tbase = mk*mk + mi*mi - 2*e_CM(s,mk,mx)*e_CM(s,mi,mj);

	// Factors (independent of cos_t) that are used for the calculation of energies Ei and Ej in the plasma frame (see below)
	double A1 = e_CM(s,mi,mj)*y;
	double A2 = e_CM(s,mj,mi)*y; //(note that the masses in e_CM have a different order)
	double A3 = pi_CM*( q_CM - q*cos_xCM )/esum;
	double A4 = pi_CM*q*sqrt(1.0 - cos_xCM*cos_xCM)/esum;


	// Add the following quantities to the array of physical parameters
	g_p.pars.insert(g_p.pars.end(), {s,q_CM,pi_CM,tbase,A1,A2,A3,A4});
	

	double integral_cos_t; 


	gsl_function F_cos_t;
	F_cos_t.function = &integrand_cos_t;
	F_cos_t.params = &g_p;

	// cos_t_min/max is always {-1,1}
	gsl_integration_qag (&F_cos_t, -1.0, 1.0, 0, rel_acc_cos_t, mem_alloc_cos_t, g_p.ipars->gsl_GK_key_cos_t, g_p.ipars->w_cos_t, &integral_cos_t, &error_bulk);

	return integral_cos_t*pi_CM/sqrt(s);

}

/*
									  _
								   | |
								  _| |_
								  \   /
								   \ /
								    '
*/


// cos_t is the (cosine of the) angle between p_x_CM and p_i_CM (in the CM frame) (t invariant angle)
double integrand_cos_t (double cos_t, void * params)
{
	
	// Create a pointer to the structure of parameters (because we don't use all of the numerical parameters and we can only dereference some of them below)
	struct general_parameters *g_p = static_cast<general_parameters *>(params);
	
	// pars = {mi,mj,mk,mx,T^2,q,ex,p,ek,sbase,esum,s,q_CM,pi_CM,tbase,A1,A2,A3,A4};

	// References to the values of parameters pointed by g_p
	// const double &mi = g_p.pars[0];
	// const double &mj = g_p.pars[1];
	// const double &mk = g_p.pars[2];
	// const double &mx = g_p.pars[3];
	const double &T2 = g_p->pars[4]; // temperature squared
	// const double &q = g_p.pars[5];
	// const double &ex = g_p.pars[6];
	// const double &p = g_p.pars[7];
	// const double &ek = g_p.pars[8];
	//const double &sbase = g_p.pars[9];
	const double &esum = g_p->pars[10];
	const double &s = g_p->pars[11];
	const double &q_CM = g_p->pars[12];
	const double &pi_CM = g_p->pars[13];
	const double &tbase = g_p->pars[14];
	const double &A1 = g_p->pars[15];
	const double &A2 = g_p->pars[16];
	const double &A3 = g_p->pars[17];
	const double &A4 = g_p->pars[18];

	// Pointer to the integration parameters
	struct integration_parameters * ipars = g_p->ipars;

	
	// Mandelstam t variable (in CM frame)
	double t = tbase + 2*q_CM*pi_CM*cos_t;

	
	//Factors for calculation of Ei and Ej in the plasma frame (independent of cos_phi) (see below)
	double E1 = A1 + A3*cos_t;
	double E2 = A2 - A3*cos_t;
	double E3 = A4*sqrt(1.0 - cos_t*cos_t);

	
	double integral_cosphi;
	
	
	// THIS IS JUST SOME TEMPORARY ADAPTATIONS (WITHOUT ANY DEFINITE REASONING) TO INCREASE THE SPEED OF THE COMPUTATION
	// CAN BE SWITCHED OFF
	if (E1 > 5 && E2 > 5 && E1 > E3)
 	{
		integral_cosphi = M_PI*exp( -esum );
	}
	else if (E1 < 5 && E2 < 5 && E1 > E3 && E3 < 1E-2) // JUST A TEMPORARY SOLUTION THAT BOOSTS THE SPEED BY A FACTOR OF 2 OR SO
	{
		integral_cosphi = M_PI/( (exp(E1)+sign_fi)*(exp(E2)+sign_fj) );
	}
	else if (E1 > 5 && E2 < 5 && E3 < 0.1) // 0.1 - ALSO A TEMPORARY SOLUTION
	{
		integral_cosphi = M_PI*exp(-E1)/sqrt(1-2*exp(E2)-(E3*E3-1)*exp(2*E2));
	}
	else // MAIN CONIDITION - EVERYTHING ABOVE CAN BE SWITCHED OFF
	{
		
		// Create a double vector of parameters for cos_phi integrand (just three values)
		std::vector<double> par_phi = {E1,E2,E3};

		gsl_function F_cos_phi;
		F_cos_phi.function = &integrand_cosphi;
		F_cos_phi.params = &par_phi;

		gsl_integration_qag (&F_cos_phi, -1.0, 1.0, 0, rel_acc_cos_phi, mem_alloc_cos_phi, ipars->gsl_GK_key_cos_phi, ipars->w_cos_phi, &integral_cosphi, &error_bulk);
	} 

	return Msquared(s*T2,t*T2)*integral_cosphi;
	
}

/*
									  _
								   | |
								  _| |_
								  \   /
								   \ /
								    '
*/

// cosphi is the (cosine of the) angle between the projections of p_x and p_i_CM onto the plane that is orthogonal to p_x_CM
double integrand_cosphi (double cosphi, void * params)
{
	std::vector<double> par_phi = * (std::vector<double> *) params;

	return 1.0/( (exp( par_phi[0] - par_phi[2]*cosphi )+sign_fi)*(exp( par_phi[1] + par_phi[2]*cosphi )+sign_fj)*sqrt(1.0 - cosphi*cosphi) );
}
// We try to make integrand_cosphi as short as possible since it is evaluated the most number of times 

