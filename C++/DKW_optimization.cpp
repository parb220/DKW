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
#include "CEquiEnergyModel_CAR.hpp"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "CMetropolis.h"
#include "TaskScheduling.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std; 

int main(int argc, char **argv)
{
	// Parallelization
	MPI_Init(&argc, &argv);
        int my_rank, nNode;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nNode);
        dw_initialize_generator(time(NULL));


	static struct option long_options[] =
	{
		{"bound file", required_argument, 0, 'b'}, 	
		// Equi-energy simulation option
		{"run id", required_argument, 0, 'r'},
                {"Number of striations", required_argument, 0, 'N'},
                {"nInitial", required_argument, 0, 'S'},
		{"percentile", required_argument, 0, 'P'},
		// Maximization option
		{"perturbation iteration", required_argument, 0, 'e'},
		{"perturbation scale", required_argument, 0, 'E'},
		{"maximization iteration", required_argument, 0, 'M'}, 
		{0, 0, 0, 0}	
	};	
	
	int option_index = 0; 
	int nInitial = 1; 
	
	CEESParameter sim_option; 
	sim_option.storage_marker = 10000;
        sim_option.run_id = to_string(time(NULL));
	sim_option.storage_dir = getenv("HOME")+string("/MATLAB_FRBA_BIN/work/results/");
	sim_option.lambda_1 = 0.1;
        sim_option.THIN = 50;
        sim_option.pee = 1.0/(10.0*sim_option.THIN);
        sim_option.simulation_length = 200000;
        sim_option.number_energy_stage = sim_option.number_striation = 1;
        sim_option.highest_stage = sim_option.lowest_stage = -1;

        double PerturbationScale = 1.0;
        int MaxPerturbationIterations = 10;
        int MaxOptimizationIterations = 10;

	string 	p_filename, bound_filename; 
	double prctl = 0.10;
	
	while(1)
	{
		int c = getopt_long(argc, argv, "b:r:N:S:P:e:E:M:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'b':
				bound_filename = string(optarg); break; 
			case 'r': 
				sim_option.run_id = string(optarg); break; 
                        case 'N':
                                sim_option.number_striation = atoi(optarg); break;
			case 'S':
				nInitial = atoi(optarg); break; 
                        case 'P':
                                prctl = atof(optarg); break;
                        case 'e':
                                PerturbationScale = atof(optarg); break;
                        case 'E':
                                MaxPerturbationIterations = atoi(optarg); break;
                        case 'M':
                                MaxOptimizationIterations = atoi(optarg); break;
			default:
				break; 
		}
	}
	sim_option.SetTemperature_geometric();

	//////////////////////////////////////////////
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

	////////////////////////////////
	CEquiEnergyModel_CAR simulation_model; 
	simulation_model.target_model = &model; 
	simulation_model.timer_when_started = -1;
	simulation_model.metropolis = new CMetropolis(&simulation_model);
        simulation_model.parameter = &sim_option;
	simulation_model.storage = new CStorageHead(my_rank, sim_option.run_id, sim_option.storage_marker, sim_option.storage_dir, sim_option.number_energy_stage);

	// Initial guess of parameters
	TDenseVector variable_p(model.NumberVariableParameters(),0.0); 
	variable_p.RandomNormal(); 
	simulation_model.current_sample = CSampleIDWeight(variable_p); 
	CSampleIDWeight mode = simulation_model.current_sample; 
	
	// bound file	
	if (!bound_filename.empty())
	{
		ifstream input_file; 
		input_file.open(bound_filename.c_str()); 
		if (!input_file)
			cerr << "Error in opening " << bound_filename << endl; 
		else 
		{
			TDenseMatrix bound(model.NumberVariableParameters(),2,0.0); 
			input_file >> bound; 
			input_file.close(); 
			simulation_model.lBound = bound.ColumnVector(0); 
			simulation_model.uBound = bound.ColumnVector(1); 
		}
	}

	vector<TIndex> blocks; 
	blocks.push_back(TIndex(0,model.NumberVariableParameters()-1)); 

	if (my_rank == 0)
	{
		if (!simulation_model.storage->makedir())
                {
                        cerr << "Error in making directory for " << sim_option.run_id << endl;
                        double *sMessage= new double [N_MESSAGE];
                        for (int i=1; i<nNode; i++)
                                MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
                        delete [] sMessage;
                        exit(1);
                }
		master_mode_finding_deploying(nNode, nInitial, simulation_model, mode, prctl);
	}
	else 
                slave_mode_finding_computing(nInitial, simulation_model, mode, MaxOptimizationIterations, MaxPerturbationIterations, PerturbationScale);
}
