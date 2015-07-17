#ifndef _CAR_HEADER_
#define _CAR_HEADER_

#include "dw_dense_matrix.hpp"

class TDenseMatrix; 
class TDenseVector; 
class CDatedDataSeries; 
class CData_FRBA; 
class CAR; 
class CAR_MinusLogLikelihood_NPSOL; 

class CAR_MinusLogLikelihood_NPSOL
{
public:
	static CAR *model; 
	static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate); 
}; 

class CAR
{
public:
	int Nfac;	// number of factors
	CData_FRBA *dataP; 	// pointer to data 

protected:
	TDenseVector A; 
	TDenseMatrix B; 
	TDenseMatrix OMEGA; 
	TDenseVector x0_0; 
	TDenseMatrix P0_0; 

	TDenseVector a; 
	TDenseMatrix b; 
	TDenseVector y_obs; 
	TDenseMatrix R; 
	
	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1)=0;  

public:
	double LogLikelihood(TDenseVector &logL, TDenseMatrix &state, TDenseMatrix &model_IE_options); 
	double LogLikelihood(TDenseVector &logL, TDenseMatrix &state, TDenseMatrix &model_IE_options, const TDenseVector &parameter); 

	virtual int NumberParameters() const = 0; 

 	virtual bool SetParameters(const TDenseVector &parameter) = 0;
        virtual TDenseVector GetParameters() = 0;
	
	virtual double MaximizeLogLikelihood(TDenseVector &x_0, const TDenseVector &lower_x, const TDenseVector &upper_x, const std::vector<std::string> &option, int max_perturbation_iteration=10, double perturbation_scale=1.0, int max_optimization=1); 

	// construction & destruction
	CAR(int _Nfac=0, CData_FRBA* _dataP=NULL); 
	CAR(const CAR &right);  
	~CAR(){}

	friend class CAR_MinusLogLikelihood_NPSOL; 
}; 

void YieldFacLoad(TDenseVector &ay, TDenseMatrix &by, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, double rho0, const TDenseVector &rho1, const TDenseVector &lambda0, const TDenseMatrix &SIGMAlambda1, const TDenseVector &MATgrid);

void YieldFacLoad_ODE(TDenseVector &A, TDenseMatrix &B, const TDenseVector &TAUgrid, const TDenseVector &k, const TDenseMatrix &K, const TDenseMatrix &H, double rho0, const TDenseVector &rho1);

void InfExpFacLoad(double &A, TDenseVector &B, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, const TDenseVector &sigq, double sigqx, double rho0, const TDenseVector &rho1, double maturity); 

TIndex FixedVariableParameter(TDenseMatrix &destination, const TDenseMatrix & source, int offset); 
TIndex FixedVariableParameter(TDenseVector &destination, const TDenseVector &source, int offset); 

bool iequals(const std::string& a, const std::string& b);
#endif
