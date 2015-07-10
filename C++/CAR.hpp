#ifndef _CAR_HEADER_
#define _CAR_HEADER_

#include "dw_dense_matrix.hpp"

class TDenseMatrix; 
class TDenseVector; 
class CDatedDataSeries; 
class CData_FRBA; 

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
	
	virtual void ObservationEquation(int j)=0;  

public:
	double LogLikelihood(TDenseVector &logL, TDenseMatrix &state); 
	double LogLikelihood(TDenseVector &logL, TDenseMatrix &state, const TDenseVector &parameter); 

	virtual int NumberParameters() const = 0; 

 	virtual bool SetParameters(const TDenseVector &parameter) = 0;
        virtual TDenseVector GetParameters() = 0;

	// construction & destruction
	CAR(int _Nfac=0, CData_FRBA* _dataP=NULL); 
	CAR(const CAR &right);  
	~CAR(){}
}; 

void YieldFacLoad(TDenseVector &ay, TDenseMatrix &by, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, double rho0, const TDenseVector &rho1, const TDenseVector &lambda0, const TDenseMatrix &SIGMAlambda1, const TDenseVector &MATgrid);

void YieldFacLoad_ODE(TDenseVector &A, TDenseMatrix &B, const TDenseVector &TAUgrid, const TDenseVector &k, const TDenseMatrix &K, const TDenseMatrix &H, double rho0, const TDenseVector &rho1);

#endif
