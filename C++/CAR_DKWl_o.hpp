#ifndef _CAR_DKWl_o_HEADER_
#define _CAR_DKWl_o_HEADER_

#include "CAR_DKWl.hpp"
#include "dw_dense_matrix.hpp"


class CAR_DKWl_o : public CAR_DKWl
{
protected:
	TDenseVector delta_options; 

	TDenseVector aI_Q; 
	TDenseMatrix bI_Q; 

	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
	virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);	
	virtual bool SetParameters_InitializeABOmega(); 
	using CAR_DKWl::UpdateABOmega;  
public:
	using CAR_DKWl:: SetAsFixed; 
	virtual bool SetAsFixed(const TDenseVector &v, const std::string &which_one);

	virtual int NumberParameters() const; 	
	virtual TDenseVector GetParameters(); 
	CAR_DKWl_o(int _Nfac=0, CData_FRBA *_dataP=NULL); 
	CAR_DKWl_o(const CAR_DKWl_o &right); 
	virtual ~CAR_DKWl_o() {}; 
};

#endif
