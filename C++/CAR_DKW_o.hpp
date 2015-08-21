#ifndef _DKW_O_HEADER_
#define _DKW_O_HEADER_

#include "CAR_DKW.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKW_o : public CAR_DKW
{
protected:
	TDenseVector delta_options; 

	TDenseVector aI_Q; 
	TDenseMatrix bI_Q; 

        virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
        virtual bool SetParameters_InitializeABOmega(); 
        using CAR_DKW::UpdateABOmega; 
public:
        using CAR_DKW::SetAsFixed;
        virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);
        virtual bool SetAsFixed(double v, const std::string &which_one);

        virtual int NumberParameters() const;
        virtual TDenseVector GetParameters();

        CAR_DKW_o(int _Nfac=0, CData_FRBA *_dataP=NULL);
        CAR_DKW_o(const CAR_DKWl &right);
        virtual ~CAR_DKW_o() {};
}; 
#endif
