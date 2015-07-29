#ifndef _CAR_DKWl_o_HEADER_
#define _CAR_DKWl_o_HEADER_

#include "CAR_DKW.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKWl_o : public CAR_DKW
{
protected:
	TDenseVector delta_options; 
	double KAPPA_L;
	double SIGMA_L; 
	double theta_L; 
	TDenseVector rho1_L; 
	double rhoL_L;
	double lambda0_L;
	double SIGMAlambda1_L;

	TDenseVector ay_TR; 
	TDenseMatrix by_TR; 
	TDenseVector ay_L; 
	TDenseMatrix by_L; 
	TDenseVector aI_Q; 
	TDenseMatrix bI_Q; 

	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
	virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);	
	virtual bool SetParameters_InitializeABOmega(); 
	virtual bool UpdateABOmega(const TDenseVector &xtm1_tm1) { return true; }; 
public:
	virtual bool SetAsFixed(const TDenseMatrix &M, const std::string &which_one)
	{
		return CAR_DKW::SetAsFixed(M,which_one);
	}
	virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);
        virtual bool SetAsFixed(double v, const std::string &which_one);

	virtual int NumberParameters() const; 	
	virtual TDenseVector GetParameters(); 
	CAR_DKWl_o(int _Nfac=0, CData_FRBA *_dataP=NULL); 
	CAR_DKWl_o(const CAR_DKWl_o &right); 
	~CAR_DKWl_o() {}; 
};

#endif
