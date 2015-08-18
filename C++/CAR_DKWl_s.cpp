#include "dw_dense_matrix.hpp"
#include "prcsn.h"
#include "matrix_operation.hpp"
#include "CAR_DKWl_s.hpp"
#include "CData_FRBA.hpp"

CAR_DKWl_s::CAR_DKWl_s(int _Nfac, CData_FRBA *_dataP) :
CAR_DKWl(_Nfac, _dataP), 
delta_swap() 
{
	if (dataP)
		delta_swap = TDenseVector(dataP->ytips_vec.cols,MINUS_INFINITY); 
}

CAR_DKWl_s::CAR_DKWl_s(const CAR_DKWl_s &right) : 
CAR_DKWl(right), 
delta_swap(right.delta_swap)
{
	delta_swap.CopyContent(right.delta_swap); 
}

TDenseVector CAR_DKWl_s::GetParameters()
{
	if (!if_internal_parameter_set)
	{
		TDenseVector complete_parameter(NumberParameters(),0.0);
		complete_parameter.Insert(0, CAR_DKWl::GetParameters());
                int counter = CAR_DKWl::NumberParameters();

		if (delta_swap.Dimension() != dataP->ytips_vec.cols)
			delta_swap = TDenseVector(dataP->ytips_vec.cols, MINUS_INFINITY); 
		complete_parameter.Insert(counter, delta_swap); 
		counter += delta_swap.Dimension(); 

		internal_parameter.CopyContent(complete_parameter); 
		if_internal_parameter_set = true; 
	}
	return internal_parameter; 
}

bool CAR_DKWl_s :: SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (_parameter.Dimension() != CAR_DKWl_s::NumberParameters())
		return false; 
	if (!CAR_DKWl::SetParameters_FromVectorToMatrix(_parameter.SubVector(0,CAR_DKWl::NumberParameters()-1)) )
		return false; 

	int counter = CAR_DKWl::NumberParameters(); 
	if (delta_swap.Dimension() != dataP->ytips_vec.cols)
		delta_swap.Resize(dataP->ytips_vec.cols); 
	delta_swap.Insert(0,_parameter, counter, counter+delta_swap.Dimension()-1); 
	counter += delta_swap.Dimension(); 

	return true; 
}

TDenseVector CAR_DKWl_s::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	TDenseVector model_IE_options_vector = CAR_DKWl::ObservationEquation(j, xt_tm1);
	for (int kk=0; kk<dataP->swap_vec.cols; kk++)
	{
		if (dataP->swap_vec(j, kk) != MINUS_INFINITY)
		{
			TDenseMatrix tmp_b(b); 
			tmp_b.CopyContent(b); 
			b.Zeros(tmp_b.rows+1, tmp_b.cols); 
			b.Insert(0,0,tmp_b); 
			b.InsertRowMatrix(tmp_b.rows,1,by_R.RowVector(kk)); 

			TDenseVector tmp_a(a); 
			tmp_a.CopyContent(a); 
			a.Zeros(tmp_a.Dimension()+1); 
			a.Insert(0,tmp_a); 
			a(tmp_a.Dimension()) = ay_R(kk); 

			TDenseVector tmp_y_obs(y_obs); 
			tmp_y_obs.CopyContent(y_obs); 
			y_obs.Zeros(tmp_y_obs.Dimension()+1); 
			y_obs.Insert(0,tmp_y_obs); 
			y_obs(tmp_y_obs.Dimension()) = dataP->yS_vec(j,kk); 

			TDenseMatrix tmp_R(R); 
			tmp_R.CopyContent(R); 
			R.Zeros(tmp_R.rows+1,tmp_R.cols+1); 
			R.Insert(0,0,tmp_R); 
			R(tmp_R.rows, tmp_R.cols) = delta_swap(kk)*delta_swap(kk); 
		}
	}
	return model_IE_options_vector; 
}

bool CAR_DKWl_s:: SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
	TIndex local_fixed_index;
        int offset;
        if (iequals(which_one, std::string("delta_swap")) )
        {
                offset = CAR_DKWl::NumberParameters();
                local_fixed_index = FixedVariableParameter(delta_swap, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else
                return CAR_DKWl::SetAsFixed(v, which_one);

}

int CAR_DKWl_s::NumberParameters() const
{
	return CAR_DKWl::NumberParameters() + dataP->ytips_vec.cols; 
}

