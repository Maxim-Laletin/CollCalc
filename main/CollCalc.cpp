#include <CollCalc.h>

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

// We use this function to fill the parameter vectors from the corresponding files
void fill_vector (std::vector<double> * vec, std::string filename)
{	
	std::ifstream infile(filename); // Open the input file

    if (infile.is_open()) 
    {
    	std::string line;
    	while (std::getline(infile, line,',')) // if there is a line
    	{ 
	        std::istringstream iss(line); // Read a line from the file

	        double num;

	        // Read numbers from the line into the vector
	        while (iss >> num) { vec->push_back(num); }
    	}
    	//else { printf("Error: Empty file.\n");  }  
    }
    else{ printf("Error: Unable to open input file.\n"); }

    infile.close();
}


// Mass of the heaviest particle in the reaction
double Mmax;
// Signs in front of the PDF (depending on the statistics)
double sign_fi, sign_fj, sign_fk;
// Desired relative accuracies for each integration
double rel_acc_ek, rel_acc_cos_s, rel_acc_cos_t, rel_acc_cos_phi; 

int mem_alloc_ek, mem_alloc_cos_s, mem_alloc_cos_t, mem_alloc_cos_phi;
int gsl_GK_key_ek, gsl_GK_key_cos_s, gsl_GK_key_cos_t, gsl_GK_key_cos_phi;


// ******************************************************************

int main(int argc, char* argv[]) 
{
	double xmin, xmax, qmin, qmax;
	int  Nx, Nq;
	//double x, q; // x is defined as x = m/T, where m is the mass of the heaviest particle in the reaction (mi), q = p/T
	

	double rel_err, result; // relative error of the integration and result

	std::string input_filepath = "parameters";
	std::string output_filepath = "bin";

	// Find the mass of the heaviest particle (to define x)
	Mmax = std::max(std::max(Mi,Mj),std::max(Mk,Mx));

	// Assign the signs to the PDFs
	sign_fi = assign_sign(s_i);
	sign_fj = assign_sign(s_j);
	sign_fk = assign_sign(s_k);

	
	// !!! CHANGE THE NAMES OF THESE PARAMETERS AND TRY TO UNIFY THEIR USAGE !!!
	
	// Desired relative accuracy of the integration
	rel_acc_cos_s = 0.1;
	rel_acc_cos_t = 0.1;
	rel_acc_cos_phi = 0.1;

	// Size of memory allocated for each integration
	mem_alloc_ek = 1000;
	mem_alloc_cos_s = 1000;
	mem_alloc_cos_t = 1000;
	mem_alloc_cos_phi = 1000;

	// Keys that correspond to different number of nodes in the GK quadrature
	gsl_GK_key_ek = 2;
	gsl_GK_key_cos_s = 1;
	gsl_GK_key_cos_t = 1;
	gsl_GK_key_cos_phi = 1;


	// The function-pointer points to the function M2(s,t) defined in the model.cpp file
	Msquared = M2;


	#ifdef INT3
	// For a I3 integral we need to provide 3 physical parameters: x,q and p + one numerical parameter that is the desired accuracy

		if (argc == 5) // 4 arguments in the command line
		{
			double x = atof(argv[1]); // [x]=1
			double q = atof(argv[2]); // [q]=1
			double p = atof(argv[3]); // [p]=1
			rel_acc_ek = atof(argv[3]); // relative accuracy

			result = integral_cos_s(x,q,p,rel_err);

			printf("result = %.3e (1.0 ± %.3e) \n",result,rel_err);

		}
		else
		{
			printf("Wrong number of input arguments!\n Should be 3 {x,q,p} \n");
		}

		// ONE CAN ADD HERE SOMETHING SIMILAR TO WHAT IS DONE IN THE INT4 PIECE OF CODE BELOW - RUNNING IN PARALLEL OVER RANGE OF X, Q AND P WHEN JUST ONE ARGUMENT IS PROVIDED


	#elif defined(INT4)
	// For a I4 integral we need to provide 2 physical parameters: x and q + one numerical parameter that is the desired accuracy


		if (argc == 4) // 3 arguments in the command line
		{
			double x = atof(argv[1]); // [x]=1
			double q = atof(argv[2]); // [q]=1
			rel_acc_ek = atof(argv[3]); // relative accuracy

			result = integral_ek(x,q,rel_err);

			printf("result = %.3e (1.0 ± %.3e) \n",result,rel_err);

		}
		else if (argc == 2) // 1 argument in the command line
		{
			rel_acc_ek = atof(argv[1]); // relative accuracy
			
			std::vector<double> x,q;

			// We take these values from the csv files in the "parameters" folder
			fill_vector(&x,input_filepath + "/x.csv");
			fill_vector(&q,input_filepath + "/q.csv");

			// Sizes of vectors
			Nx = x.size();
			Nq = q.size();

			// Declare a 2D array
			std::vector<std::vector<double>> result(Nx, std::vector<double>(Nq));

			// Running the computation over the range of x and q values in parallel (!) 
			#pragma omp parallel for schedule(dynamic,1) collapse(2)
			for (int i = 0; i < Nx; i++) {
			    for (int j = 0; j < Nq; j++) {
			        result[i][j] = integral_ek(x[i],q[j],rel_err);
			    }
			}

			// Writing to file
			char filename[64];
			snprintf(filename, sizeof(filename), "bin/CollTerm%s.csv", pr_name.c_str());


		    FILE *fp; // initiate a file
		    fp = fopen(filename,"w"); // open the file

		    // Writing the table to the file using nested loops
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Nq; j++)
				{	
					fprintf(fp,"%.3e ",result[i][j]);
				}
				fprintf(fp,"\n");
			}

			fclose(fp);


		}
		else 
		{
			printf("Wrong number of input arguments!\n Should be either 2 {x,q} or 6 {xmin,xmax,qmin,qmax,Nx,Nq}\n");
		}

	#endif

	return 0;

}
