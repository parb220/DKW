#include <vector>
#include <ctime>
#include "dw_dense_matrix.hpp"
#include "dw_rand.h"
#include "prcsn.h"
#include "CEquiEnergyModel_CAR.hpp"
#include "CAR.hpp"

using namespace std; 

double CEquiEnergyModel_CAR::log_posterior_function(CSampleIDWeight &x)
{
	if (!x.calculated)
	{
		TDenseVector logL; 
		TDenseMatrix state, model_IE_options; 
		x.reserved = target_model->LogLikelihood(logL, state, model_IE_options, x.data);
		x.weight = x.reserved + log_prior_function(x); 
		x.calculated = true; 
	}
	if (if_bounded)
		return x.weight * lambda; 
	else 
		return x.weight; 
}

double CEquiEnergyModel_CAR::log_likelihood_function(const CSampleIDWeight &x)
{
	TDenseVector logL;
 	TDenseMatrix state, model_IE_options;
	return target_model->LogLikelihood(logL, state, model_IE_options, x.data); 
}

double CEquiEnergyModel_CAR::log_prior_function(const CSampleIDWeight &x)
{
	if (lBound.Dimension() == 0 && uBound.Dimension() == 0)
		return 0;  
	else if (lBound.Dimension() == x.data.Dimension() && uBound.Dimension() == 0)
	{
		for (int i=0; i<lBound.Dimension(); i++)
		{
			if (x.data(i) < lBound(i)) 
				return MINUS_INFINITY; 
		}
		return 0; 
	}
	else if (lBound.Dimension() == 0 && uBound.Dimension() == x.data.Dimension())
	{
		for (int i=0; i<uBound.Dimension(); i++)
		{
			if (x.data(i) > uBound(i))
				return MINUS_INFINITY; 
		}
		return 0; 
	}
	else if (lBound.Dimension() == x.data.Dimension() && uBound.Dimension() == x.data.Dimension())
	{
		for (int i=0; i<x.data.Dimension(); i++)
		{
			if (x.data(i) < lBound(i) || x.data(i) > uBound(i))
				return MINUS_INFINITY; 
		}
		return 0; 
	}
	else 
		return MINUS_INFINITY; 
}

double CEquiEnergyModel_CAR::log_posterior_function(const double *x, int n)
{
	TDenseVector x_vec(n); 
	memcpy(x_vec.vector, x, sizeof(double)*n); 
	CSampleIDWeight sample(x_vec); 
	return log_posterior_function(sample); 	
}

double CEquiEnergyModel_CAR::log_likelihood_function(const double *x, int n)
{
	TDenseVector x_vec(n); 
	memcpy(x_vec.vector, x, sizeof(double)*n); 
	CSampleIDWeight sample(x_vec); 
	return log_likelihood_function(sample); 
}

double CEquiEnergyModel_CAR::log_prior_function(const double *x, int n)
{
	TDenseVector x_vec(n);
        memcpy(x_vec.vector, x, sizeof(double)*n);
	CSampleIDWeight sample(x_vec); 
        return log_prior_function(sample);
}

bool CEquiEnergyModel_CAR::DrawParametersFromPrior(double *x) const
{
	TDenseVector x_vec(current_sample.data.Dimension()); 
	x_vec.RandomNormal(); 
	memcpy(x,x_vec.vector, sizeof(double)*x_vec.Dimension()); 
	return true; 		
}

CEquiEnergyModel_CAR::CEquiEnergyModel_CAR() :
CEquiEnergyModel(), 
target_model(NULL), 
lBound(),
uBound()
{}

CEquiEnergyModel_CAR::CEquiEnergyModel_CAR(bool _if_bounded, int eL, double _lambda, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage, CAR *_model, const TDenseVector &_lBound, const TDenseVector &_uBound) :
CEquiEnergyModel(_if_bounded, eL, _lambda, _x, _time, _metropolis, _parameter, _storage),
target_model(_model), lBound(_lBound), uBound(_uBound)
{
	lBound.CopyContent(_lBound); 
	uBound.CopyContent(_uBound); 
}

CEquiEnergyModel_CAR::~CEquiEnergyModel_CAR()
{}
