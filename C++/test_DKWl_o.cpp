#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <getopt.h>
#include "prcsn.h"
#include "CDate.hpp"
#include "CData_FRBA.hpp"
#include "CAR_DKWl_o.hpp"

using namespace std; 

int main(int argc, char **argv)
{

	static struct option long_options[] =
        {
                {"initial parameter file", required_argument, 0, 'p'},
                {"bound file", required_argument, 0, 'b'},
                {"perturbation iteration", required_argument, 0, 't'},
                {"perturbation scale", required_argument, 0, 's'},
                {"maximization iteration", required_argument, 0, 'm'},
                {0, 0, 0, 0}
        };

	string  p_filename, bound_filename;
        int perturbation_iteration = 10, maximization_iteration = 1;
        double perturbation_scale = 1.0;

	int option_index = 0;

        while(1)
        {
                int c = getopt_long(argc, argv, "p:b:t:s:m:", long_options, &option_index);
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
                        default:
                                break;
                }
        }

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

	int Nfac = 3; 
	CAR_DKWl_o model(Nfac, &data); 
	TDenseMatrix _kappa = DiagonalMatrix(MINUS_INFINITY, Nfac, Nfac); 
	model.SetAsFixed(_kappa, string("KAPPA")); 
	
	TDenseMatrix _sigma(Nfac, Nfac, MINUS_INFINITY); 
	_sigma(0,0) = 0.01; _sigma(0,1) = 0.0; _sigma(0,2) = 0.0; 
	_sigma(1,1) = 0.01; _sigma(1,2) = 0.0;
	_sigma(2,2) = 0.01; 
	model.SetAsFixed(_sigma, string("SIGMA")); 

	model.SetAsFixed(0.01, string("SIGMA_L")); 

	TDenseVector _theta(Nfac,0.0); 
	model.SetAsFixed(_theta, string("theta")); 

	model.SetAsFixed(0.0075, string("delta_bcfLT")); 
	model.SetAsFixed(0.0, string("delta_p")); 

	TDenseVector variable_p(model.NumberVariableParameters(),0.0); 
	variable_p.RandomNormal(); 
	if (!p_filename.empty())
	{
		ifstream input_file; 
		input_file.open(p_filename.c_str()); 
		if (!input_file)
		{
			cerr << "Error in opening " << p_filename << endl; 
			exit(-1);
		}
		input_file >> variable_p; 
		input_file.close(); 
	}
	
	TDenseMatrix bound(model.NumberVariableParameters(),2,0.0); 
	bound.InsertColumnMatrix(0, 0, MINUS_INFINITY*Zeros(model.NumberVariableParameters())); 
	bound.InsertColumnMatrix(0, 1, PLUS_INFINITY*Zeros(model.NumberVariableParameters())); 
	if (!bound_filename.empty())
	{
		ifstream input_file; 
		input_file.open(bound_filename.c_str()); 
		if (!input_file)
		{
			cerr << "Error in opening " << bound_filename << endl; 
			exit(-1); 
		}
		input_file >> bound; 
		input_file.close();
	}

	TDenseVector logLSeries; 
	TDenseMatrix state;
	TDenseMatrix model_IE_options;  
	double log_likelihood = model.LogLikelihood(logLSeries, state, model_IE_options, variable_p); 
	cout << setprecision(20) << log_likelihood; 
	for (int i=0; i<variable_p.Dimension(); i++)
		cout << "\t" << variable_p[i]; 
	cout << endl; 
	
	double optimal_log_likelihood = model.MaximizeLogLikelihood(variable_p, TDenseVector(), TDenseVector(), std::vector<std::string>(0), perturbation_iteration, perturbation_scale, maximization_iteration); 
	cout << setprecision(20) << optimal_log_likelihood; 
	for (int i=0; i<variable_p.Dimension(); i++)
		cout << "\t" << variable_p[i]; 
	cout << endl; 
}
