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
delta_options(right.delta_options), aI_Q(), bI_Q()
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
	if (!CAR_DKW::SetParameters_FromVectorToMatrix(_parameter.SubVector(0,CAR_DKW::NumberParameters()-1))
		return false; 

	int counter = CAR_DKW::NumberParameters(); 
	if (delta_options.Dimension() != dataP->MATgrid_options.Dimension())
		delta_options.Resize(dataP->MATgrid_options.Dimension()); 
	delta_options.Insert(0, parameter, counter, counter+delta_options.Dimension()-1); 
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
                TDenseVector temp_by_vector = -MAT * temp_by.RowVector(0);

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
	TDenseVector model_IE_options = CAR_DKW::ObservationEquation(j, xt_tm1); 
	if (dataP->IE_options(j,0) != MINUS_INFINITY)
	{
		
	}
}
