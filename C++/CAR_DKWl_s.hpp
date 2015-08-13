#ifndef _CAR_DKWL_S_HEADER
#define _CAR_DKWL_S_HEADER

#include "CAR_DKWl.hpp"
#include "dw_dense_matrix.hpp"

class CAR_DKWl_s : public CAR_DKWl
{
protected: 
	TDenseVector delta_swap; 

	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
        virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);
	using CAR_DKWl::SetParameters_InitializeABOmega; 
	using CAR_DKWl::UpdateABOmega; 

public: 
	using CAR_DKWl:: SetAsFixed;
        virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);

        virtual int NumberParameters() const;
        virtual TDenseVector GetParameters();

	CAR_DKWl_s(int _Nfac=0, CData_FRBA *_dataP=NULL); 
	CAR_DKWl_s(const CAR_DKWl_s &right); 
	virtual ~CAR_DKWl_s() {}; 
}; 

#endif
