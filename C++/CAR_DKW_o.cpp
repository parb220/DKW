#include "CAR_DKW_o.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"
#include "prcsn.h"

CAR_DKW_o::CAR_DKW_o(int _Nfac, CData_FRBA *_dataP) :
CAR_DKW(_Nfac, _dataP), 
delta_options(), aI_Q(), bI_Q()
{
	if (dataP)
		delta_options = TDenseVector(dataP->MATgrid_options.Dimension(), MINUS_INFINITY); 
}

CAR_DKW_o::CAR_DKW_o(const CAR_DKW_o &right) :
CAR_DKW(right), 
delta_options(right.delta_options), aI_Q(right.aI_Q), bI_Q(right.bI_Q)
{
	delta_options.CopyContent(right.delta_options); 
	aI_Q.CopyContent(right.aI_Q); 
	bI_Q.CopyContent(right.bI_Q); 
}

TDenseVector CAR_DKW_o::GetParameters()
{
	if (!if_internal_parameter_set)
	{
		TDenseVector complete_parameter(NumberParameters(),0.0); 
		complete_parameter.Insert(0,CAR_DKW::GetParameters()); 
		int counter = CAR_DKW::NumberParameters(); 
		
		if (delta_options.Dimension() != dataP->MATgrid_options.Dimension() )
			delta_options = TDenseVector(dataP->MATgrid_options.Dimension(), MINUS_INFINITY); 
		complete_parameter.Insert(counter, delta_options); 
		counter += delta_options.Dimension(); 

		internal_parameter.CopyContent(complete_parameter); 
		if_internal_parameter_set = true;
	}
	return internal_parameter;
}

bool CAR_DKW_o::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (_parameter.Dimension() < NumberParameters())
		return false; 
	if (!CAR_DKW::SetParameters_FromVectorToMatrix(_parameter.SubVector(0,CAR_DKW::NumberParameters()-1)) )
		return false; 

	int counter = CAR_DKW::NumberParameters(); 
	if (delta_options.Dimension() != dataP->MATgrid_options.Dimension())
		delta_options.Resize(dataP->MATgrid_options.Dimension()); 
	delta_options.Insert(0, _parameter, counter, counter+delta_options.Dimension()-1); 
	counter += delta_options.Dimension(); 

	return true; 
}

bool CAR_DKW_o::SetParameters_InitializeABOmega()
{
	if (!CAR_DKW::SetParameters_InitializeABOmega())
		return false; 

	TDenseMatrix lambda1 = SIGMA_inverse * SIGMAlambda1;
        TDenseVector KAPPAtheta = KAPPA * theta;
        TDenseVector theta_Q;

	aI_Q.Zeros(dataP->MATgrid_options.Dimension());
        bI_Q.Zeros(dataP->MATgrid_options.Dimension(), Nfac);
        
        for (int i=0; i<dataP->MATgrid_options.Dimension(); i++)
        {
                double MAT = dataP->MATgrid_options(i);
                TDenseVector temp_ay;
                TDenseMatrix temp_by;
                if (!YieldFacLoad(temp_ay, temp_by, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, KAPPAtheta, rho0, rho1, lambda0,TDenseVector(1,MAT)))
                        return false;
                TDenseVector temp_by_vector = temp_by.RowVector(0); // -MAT * temp_by.RowVector(0);

                theta_Q = Multiply(Inv_KAPPA_rn, KAPPAtheta-SIGMA*lambda0+MultiplyTranspose(SIGMA,SIGMA)*temp_by_vector);

                double rho0_Q = rho0_pi - InnerProduct(lambda0, sigq)+InnerProduct(sigq, TransposeMultiply(SIGMA,temp_by_vector));
                TDenseVector rho1_Q = rho1_pi - TransposeMultiply(lambda1, sigq);

                double temp_aI_Q;
                TDenseVector temp_bI_Q;
                InfExpFacLoad(temp_aI_Q, temp_bI_Q, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, theta_Q, sigq, sigqx, rho0_Q, rho1_Q, MAT);

                aI_Q(i) = temp_aI_Q;
                bI_Q.InsertRowMatrix(i, 0, temp_bI_Q);
        }
	return true; 
}

TDenseVector CAR_DKW_o::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	TDenseVector model_IE_options_vector = CAR_DKW::ObservationEquation(j, xt_tm1); 
	if (dataP->IE_options(j,0) != MINUS_INFINITY)
	{
		TDenseMatrix tmp_b(b);
                tmp_b.CopyContent(b);
                b.Zeros(tmp_b.rows+bI_Q.rows, tmp_b.cols);
                b.Insert(0,0,tmp_b);
                b.Insert(tmp_b.rows, 1, bI_Q);

                TDenseVector tmp_a(a);
                tmp_a.CopyContent(a);
                a.Zeros(tmp_a.Dimension()+aI_Q.Dimension());
                a.Insert(0, tmp_a);
                a.Insert(tmp_a.Dimension(), aI_Q);

                TDenseVector tmp_y_obs(y_obs);
                tmp_y_obs.CopyContent(y_obs);
		y_obs.Zeros(tmp_y_obs.Dimension()+2);
		y_obs.Insert(0,tmp_y_obs); 
		y_obs(tmp_y_obs.Dimension()) = (dataP->IE_options(j,0)+dataP->IE_options(j,1)+dataP->IE_options(j,2))/3.0;
                y_obs(tmp_y_obs.Dimension()+1) = (dataP->IE_options(j,3)+dataP->IE_options(j,5)+dataP->IE_options(j,5))/3.0;

		TDenseMatrix tmp_R(R);
                tmp_R.CopyContent(R);
                R.Zeros(tmp_R.rows+delta_options.Dimension(), tmp_R.cols+delta_options.Dimension());
                R.Insert(0,0,tmp_R);
                for (int i=0; i<delta_options.Dimension(); i++)
                        R(i+tmp_R.rows, i+tmp_R.cols) = delta_options(i) * delta_options(i);

                model_IE_options_vector = aI_Q + bI_Q * xt_tm1.SubVector(1,3);
	}
	return model_IE_options_vector; 
}

bool CAR_DKW_o::SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
        TIndex local_fixed_index;
        int offset;
        if (iequals(which_one, std::string("delta_options")) )
        {
                offset = CAR_DKW::NumberParameters();
                local_fixed_index = FixedVariableParameter(delta_options, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else
                return CAR_DKW::SetAsFixed(v, which_one);
}

int CAR_DKW_o::NumberParameters() const
{
        return CAR_DKW::NumberParameters() + dataP->MATgrid_options.Dimension();
}


