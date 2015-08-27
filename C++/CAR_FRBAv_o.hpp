#ifndef _HEADER_CAR_FRBAV_O_
#define _HEADER_CAR_FRBAV_O_

#include "CAR_DKWv_o.hpp"
#include "dw_dense_matrix.hpp"

class CAR_FRBAv_o : public CAR_DKWv_o
{
private:
	bool delta_option_position_set; 
	int delta_option_index_start; 
	int delta_option_index_end; 

	bool SetDeltaOptionPosition(); 
// delta_options no long a TDenseVector, but a double
protected:
	virtual bool SetParameters_FromVectorToMatrix(const TDenseVector &);	
	virtual TDenseVector ObservationEquation(int j, const TDenseVector &xt_tm1);
	bool CapsFacLoad_FRBAv_o(TDenseVector &, TDenseMatrix &, TDenseVector &, TDenseVector &, const TDenseVector &xt_tm1, double vt_tm1, const TDenseVector &Maturity, const TDenseVector &Strike, bool);

public:
	virtual int NumberParameters() const; 
	virtual TDenseVector GetParameters(); 
	using CAR_DKWv_o :: SetAsFixed; 
	virtual bool SetAsFixed(double , const std::string &which_one); 

	CAR_FRBAv_o(int _Nfac=0, CData_FRBA *_dataP=NULL);
	CAR_FRBAv_o(const CAR_FRBAv_o &right); 
	virtual ~CAR_FRBAv_o() {}; 
}; 

#endif
