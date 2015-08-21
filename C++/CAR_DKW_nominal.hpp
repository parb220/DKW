#ifndef _DKW_NOMINAL_HEADER_
#define _DKW_NOMINAL_HEADER_

#include "CAR_DKW.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKW_nominal : public CAR_DKW
{
private: 
	bool position_removed_from_DKW_set; 
	TIndex position_removed_from_DKW;
	bool SetParametersNotConsidered(); 
protected:
	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
        virtual bool SetParameters_InitializeABOmega();
        using CAR_DKW::UpdateABOmega;

public:
        using CAR_DKW:: SetAsFixed;
        virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);
        virtual bool SetAsFixed(double v, const std::string &which_one);

        virtual int NumberParameters() const;
        virtual TDenseVector GetParameters();

	CAR_DKW_nominal(int tNfac=0, CData_FRBA *tDataP=NULL); 
	CAR_DKW_nominal(const CAR_DKW_nominal &right); 
	virtual ~CAR_DKW_nominal() {}; 
}; 

#endif
