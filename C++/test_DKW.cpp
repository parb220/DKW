#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
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
	ifstream input_file; 
	/*string p_filename("../C++/parameter.txt"); 
	input_file.open(p_filename.c_str()); 
	if (!input_file)
	{
		cerr << "Error in opening " << p_filename << endl; 
		exit(-1);
	}
	input_file >> variable_p; 
	input_file.close(); */
	variable_p.RandomNormal(); 
	
	TDenseMatrix bound(model.NumberVariableParameters(),2,0.0); 
	string bound_filename("../C++/bound.txt"); 
	input_file.open(bound_filename.c_str()); 
	if (!input_file)
	{
		cerr << "Error in opening " << bound_filename << endl; 
		exit(-1); 
	}
	input_file >> bound; 
	input_file.close(); 

	TDenseVector logLSeries; 
	TDenseMatrix state;
	TDenseMatrix model_IE_options;  
	double log_likelihood = model.LogLikelihood(logLSeries, state, model_IE_options, variable_p); 
	cout << setprecision(20) << log_likelihood; 
	for (int i=0; i<variable_p.Dimension(); i++)
		cout << "\t" << variable_p[i]; 
	cout << endl; 
	
	double optimal_log_likelihood = model.MaximizeLogLikelihood(variable_p, bound.ColumnVector(0), bound.ColumnVector(1), std::vector<std::string>(0), 30, 1.0, 10); 
	cout << setprecision(20) << optimal_log_likelihood; 
	for (int i=0; i<variable_p.Dimension(); i++)
		cout << "\t" << variable_p[i]; 
	cout << endl; 
}
