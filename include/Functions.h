#ifndef __FUNCTIONS__
#define __FUNCTIONS__

// *************************************** General particle physics functions ************************************************

// Energy in the CM frame
inline double e_CM (double s, double m1, double m2)
{
	return 0.5*(s + pow(m1,2) - pow(m2,2))/sqrt(s);
}

// Momentum in the CM frame
inline double p_CM (double s, double m1, double m2)
{
	return 0.5*sqrt( (s - pow(m1 + m2,2))*(s - pow(m1 - m2,2))/s );
}

// ***************************************************************************************************************************

#endif
