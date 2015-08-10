#ifndef _CEQUIENERGYMODEL_CAR_HEADER_
#define _CEQUIENERGYMODEL_CAR_HEADER_

#include "CEquiEnergyModel.h"

class CSampleIDWeight; 
class CAR; 

class CEquiEnergyModel_CAR : public CEquiEnergyModel
{
public: 
	CAR *target_model; 
	TDenseVector lBound; 
	TDenseVector uBound; 
	
	virtual double log_posterior_function(CSampleIDWeight &x);
        virtual double log_likelihood_function(const CSampleIDWeight &x);
        virtual double log_prior_function(const CSampleIDWeight &x);
        virtual double log_posterior_function(const double *x, int n);
        virtual double log_likelihood_function(const double *x, int n);
        virtual double log_prior_function(const double *x, int n);

	virtual bool DrawParametersFromPrior(double *x) const;
	using CEquiEnergyModel::HillClimb_NPSOL;

	CEquiEnergyModel_CAR(); 
	CEquiEnergyModel_CAR(bool _if_bounded, int eL, double _lambda, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage, CAR *_model, const TDenseVector &_lBound, const TDenseVector &_uBound); 
	virtual ~CEquiEnergyModel_CAR(); 	
}; 

#endif
