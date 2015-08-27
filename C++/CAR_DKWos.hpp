#ifndef _CAR_DKWos_HEADER
#define _CAR_DKWos_HEADER
#include "CAR_DKWl_s.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKWos : public CAR_DKWl_s
{
protected:
	TDenseVector delta_options; 
	
	TDenseVector aI_Q;
        TDenseMatrix bI_Q;

        virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
        virtual bool SetParameters_InitializeABOmega();
        using CAR_DKWl_s::UpdateABOmega;
public:
        using CAR_DKWl_s::SetAsFixed;
        virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);

        virtual int NumberParameters() const;
        virtual TDenseVector GetParameters();

        CAR_DKWos(int _Nfac=0, CData_FRBA *_dataP=NULL);
        CAR_DKWos(const CAR_DKWos &right);
        virtual ~CAR_DKWos() {};


}; 
#endif
