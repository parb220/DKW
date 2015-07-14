#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <mpi.h>
#include <getopt.h>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <errno.h>

#include "dw_rand.h"
#include "dw_dense_matrix.hpp"
#include "prcsn.h"
#include "CDate.hpp"
#include "CData_FRBA.hpp"
#include "CAR_DKW.hpp"

using namespace std; 

void master_deploying(const string &write_directory, int nNode, int nGroup, bool if_initial_parameter, const TDenseVector &initial_p);
void slave_computing(CAR &model, const string &write_directory, const TDenseMatrix &bound, int perturbation_iteration, double perturbation_scale, int maximization_iteration);

int main(int argc, char **argv)
{
	// Parallelization
	MPI_Init(&argc, &argv);
        int my_rank, nNode;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nNode);
        dw_initialize_generator(time(NULL));

	// program setting
	string 	p_filename, bound_filename; 
	int perturbation_iteration = 10, maximization_iteration = 1, nGroup = nNode-1; 
	double perturbation_scale = 1.0; 
	string dir = getenv("HOME")+string("/MATLAB_FRBA_BIN/work/results/"); 
	string run_id = to_string(time(NULL)); 

	static struct option long_options[] =
	{
		{"initial parameter file", required_argument, 0, 'p'},
		{"bound file", required_argument, 0, 'b'}, 	
		{"perturbation iteration", required_argument, 0, 't'},
		{"perturbation scale", required_argument, 0, 's'},
		{"maximization iteration", required_argument, 0, 'm'}, 
		{"parallelization threads", required_argument, 0, 'G'},
		{"run id", required_argument, 0, 'r'},
		{0, 0, 0, 0}	
	};	
	
	int option_index = 0; 
	
	while(1)
	{
		int c = getopt_long(argc, argv, "p:b:t:s:m:G:r:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'p':
				p_filename = string(optarg); break; 
			case 'b':
				bound_filename = string(optarg); break; 
			case 't':
				perturbation_iteration = atoi(optarg); break; 
			case 's':
				perturbation_scale = atof(optarg); break; 
			case 'm':
				maximization_iteration = atoi(optarg); break; 
			case 'G':
				nGroup = atoi(optarg); break; 
			case 'r': 
				run_id = string(optarg); break; 
			default:
				break; 
		}
	}
	string write_directory = dir + run_id;  

	// Load data file
	vector<string> data_file(13); 
	data_file[0] = string("../data/FMAyield.txt"); 
	data_file[1] = string("../data/FMAtips.txt"); 
	data_file[2] = string("../data/FMAcpi.txt"); 
	data_file[3] = string("../data/caps.txt"); 
	data_file[4] = string("../data/floors.txt"); 
	data_file[5] = string("../data/allyield.txt"); 
	data_file[6] = string("../data/oil.txt"); 
	data_file[7] = string("../data/FMAbcf.txt"); 
	data_file[8] = string("../data/FMAbcfLT.txt"); 
	data_file[9] = string("../data/FMAieS.txt"); 
	data_file[10] = string("../data/FMAieL.txt"); 
	data_file[11] = string("../data/alltips.txt"); 
	data_file[12] = string("../data/swap.txt"); 

	double Dt = 7.0/365.0; 
	CData_FRBA data(Dt); 
	CDate start_date(1990,1,1), end_date(2012,12,31); 
	if (!data.LoadDate_FRBA(data_file[0], start_date, end_date))
	{
		cerr << "Error in loading date" << endl; 
		exit(-1); 
	}
	if (!data.LoadDataSeries_FRBA(data_file, start_date, end_date))
	{
		cerr << "Error in loading data series " << endl; 
		exit(-1); 
	}

	// Model construction
	int Nfac = 3; 
	CAR_DKW model(Nfac, &data); 
	TDenseMatrix _kappa = DiagonalMatrix(MINUS_INFINITY, Nfac, Nfac); 
	model.SetAsFixed(_kappa, string("KAPPA")); 
	TDenseMatrix _sigma(Nfac, Nfac, MINUS_INFINITY); 
	_sigma(0,0) = 0.01; _sigma(0,1) = 0.0; _sigma(0,2) = 0.0; 
	_sigma(1,1) = 0.01; _sigma(1,2) = 0.0;
	_sigma(2,2) = 0.01; 
	model.SetAsFixed(_sigma, string("SIGMA")); 
	TDenseVector _theta(Nfac,0.0); 
	model.SetAsFixed(_theta, string("theta")); 
	model.SetAsFixed(0.0075, string("delta_bcfLT")); 
	model.SetAsFixed(0.0, string("delta_p")); 

	// Initial guess of parameters
	TDenseVector variable_p(model.NumberVariableParameters(),0.0); 
	variable_p.RandomNormal(); 
	bool if_initial_parameter = false; 
	if (!p_filename.empty())
	{
		ifstream input_file; 
		input_file.open(p_filename.c_str()); 
		if (!input_file)
			cerr << "Error in opening " << p_filename << endl; 
		else 
		{
			input_file >> variable_p; 
			input_file.close(); 
			if_initial_parameter = true; 
		}
	}
	
	TDenseMatrix bound(model.NumberVariableParameters(),2,0.0); 
	bound.InsertColumnMatrix(0, 0, MINUS_INFINITY*Ones(model.NumberVariableParameters())); 
	bound.InsertColumnMatrix(0, 1, PLUS_INFINITY*Ones(model.NumberVariableParameters())); 
	if (!bound_filename.empty())
	{
		ifstream input_file; 
		input_file.open(bound_filename.c_str()); 
		if (!input_file)
			cerr << "Error in opening " << bound_filename << endl; 
		else 
		{
			input_file >> bound; 
			input_file.close(); 
		}
	}

	if (my_rank == 0)
	{
		// making directory
		errno = 0; // clear error number
        	int status = mkdir(write_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        	if (status !=0 && errno != EEXIST) 
		{ 	
			cerr << "error in creating directory " << write_directory << endl; 
			exit(-1); 
		}
		master_deploying(write_directory, nNode, nGroup, if_initial_parameter, variable_p);
	}
	else 
		slave_computing(model, write_directory, bound, perturbation_iteration, perturbation_scale, maximization_iteration); 
}
