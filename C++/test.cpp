#include <string>
#include <vector>
#include <iostream>
#include "prcsn.h"
#include "CDate.hpp"
#include "CData_FRBA.hpp"
#include "CAR_DKW.hpp"

using namespace std; 

int main(int argc, char **argv)
{
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

	TDenseVector variable_p(model.NumberVariableParameters(),0.0); 
	string p_filename("../C++/parameter.txt"); 
	ifstream input_file(p_filename.c_str()); 
	if (!input_file)
	{
		cerr << "Error in opening " << p_filename << endl; 
		exit(-1);
	}
	input_file >> variable_p; 
	input_file.close(); 

	if(!model.SetParameters(variable_p))
	{
		cerr << "Error in setting parameters for the model.\n"; 
		exit(-1); 
	}
	TDenseVector logLSeries; 
	TDenseMatrix state; 
	double log_likelihood = model.LogLikelihood(logLSeries, state); 
	TDenseVector current_p = model.GetParameters(); 
}
