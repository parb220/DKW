#ifndef _HEADER_CAR_DKWLV_O_
#define _HEADER_CAR_DKWLV_O_

#include "CAR_DKW.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKWlv_o : public CAR_DKW // minus sigqx
{
private: 
 	bool position_removed_from_DKW_set; 	
	TIndex position_removed_from_DKW; 
	bool SetParametersNotConsidered(); 
protected:
	double rhov; 
	double rhov_pi; 
	
	double KAPPAv;
	double SIGMAv; 
	double thetav;
	double gamv;
	double gamvx;
	double rho;
	
	TDenseVector delta_options; 
	double KAPPA_L;
	double SIGMA_L;
	double theta_L;
	TDenseVector rho1_L;
	double rhoL_L;
	double lambda0_L;
	double SIGMAlambda1_L;

	TDenseVector cy; 
	TDenseVector cy_R; 
	TDenseVector ay_TR; 
	TDenseMatrix by_TR; 
	TDenseVector cy_TR; 
	TDenseVector ay_L; 
	TDenseMatrix by_L; 
	TDenseVector cf; 
	TDenseVector aI_Q; 
	TDenseMatrix bI_Q; 
	TDenseVector cI_Q; 
	
	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
        virtual bool SetParameters_InitializeABOmega();
	virtual bool UpdateABOmega(const TDenseVector &xtm1_tm1);

public: 
	virtual bool SetAsFixed(const TDenseMatrix &M, const std::string &which_one)
        {
                return CAR_DKW::SetAsFixed(M,which_one);
        }
        virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);
        virtual bool SetAsFixed(double v, const std::string &which_one);
	
	virtual int NumberParameters() const;
	virtual TDenseVector GetParameters();

	CAR_DKWlv_o(int _Nfac=0, CData_FRBA *_dataP=NULL); 
	CAR_DKWlv_o(const CAR_DKWlv_o &right);
	~CAR_DKWlv_o() {}; 
};
#endif
