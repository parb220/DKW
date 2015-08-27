#ifndef _HEADER_CAR_DKWV_O_
#define _HEADER_CAR_DKWV_O_

#include "CAR_DKW_o.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKWv_o : public CAR_DKW_o	// minus sigqx
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

        TDenseVector cy;
        TDenseVector cy_R;
        TDenseVector cy_TR;
        TDenseVector cf;
        TDenseVector cI_Q;

	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
        virtual bool SetParameters_InitializeABOmega();
        virtual bool UpdateABOmega(const TDenseVector &xtm1_tm1);

public:
	using CAR_DKW_o::SetAsFixed; 
	virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);
        virtual bool SetAsFixed(double v, const std::string &which_one);

	virtual int NumberParameters() const;
        virtual TDenseVector GetParameters();
	
	CAR_DKWv_o(int _Nfac=0, CData_FRBA *_dataP=NULL);
        CAR_DKWv_o(const CAR_DKWv_o &right);
        virtual ~CAR_DKWv_o() {};
}; 

#endif
