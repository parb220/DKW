#ifndef _CAR_DKW_HEADER
#define _CAR_DKW_HEADER

#include "CAR.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKW : public CAR
{
protected:
	TDenseMatrix KAPPA; 
	TDenseMatrix SIGMA; 
	TDenseVector theta;
	double rho0;
	TDenseVector rho1;
	TDenseVector lambda0;
	TDenseMatrix SIGMAlambda1; 
	TDenseVector delta_y;
	TDenseVector delta_bcf;
	double delta_bcfLT;

	double rho0_pi;
	TDenseVector rho1_pi;
	TDenseVector sigq;
	double sigqx;

	double delta_p;
	TDenseVector delta_tips;

	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1); 
	
	TDenseVector ay; 
	TDenseMatrix by; 
	TDenseVector ay_R; 
	TDenseMatrix by_R; 
	TDenseVector af; 
	TDenseMatrix bf; 
	TIndex FixedIndex;
	TIndex VariableIndex; 
	bool if_variable_index_set;  
	TDenseVector internal_parameter;
	bool if_internal_parameter_set;  	

	TDenseMatrix KAPPA_inverse; 
	TDenseMatrix SIGMA_inverse; 
	TDenseMatrix KAPPA_rn; 
	TDenseMatrix Inv_KAPPA_rn; 
	TDenseMatrix Inv_Kron_KAPPA_rn; 

	virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &); 
	virtual bool SetParameters_InitializeABOmega(); 
	virtual bool UpdateABOmega(const TDenseVector &xtm1_tm1) { return true; }; 
	TIndex GetVariableIndex(); 

public:
	virtual int NumberVariableParameters() const; 
	virtual int NumberFixedParameters() const; 

	virtual bool SetAsFixed(const TDenseMatrix &M, const std::string &which_one); 
	virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one); 
	virtual bool SetAsFixed(double v, const std::string &which_one); 

	virtual int NumberParameters() const; 
	virtual bool SetParameters(const TDenseVector &parameter); 
	virtual TDenseVector GetParameters(); 
	
	CAR_DKW(int _Nfac=0, CData_FRBA* _dataP=NULL); 
	CAR_DKW(const CAR_DKW &right); 
	virtual ~CAR_DKW() {}; 
};

#endif
