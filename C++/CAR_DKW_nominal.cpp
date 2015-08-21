#include "CAR_DKW_nominal.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"
#include "prcsn.h"

CAR_DKW_nominal::CAR_DKW_nominal(int tNfac, CData_FRBA *tDataP) : 
CAR_DKW(tNfac, tDataP), 
position_removed_from_DKW_set(false),
position_removed_from_DKW()
{
	if (dataP)
	{
		int offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols; 
		position_removed_from_DKW = TIndex(offset, offset+dataP->bcf_vec.cols); // delta_bcf, delta_bcfLT
		offset = CAR_DKW::NumberParameters()-1; 
		position_removed_from_DKW += TIndex(offset-dataP->ytips_vec.cols+1, offset); 
		position_removed_from_DKW_set = true; 
	}
}

CAR_DKW_nominal::CAR_DKW_nominal(const CAR_DKW_nominal &right) :
CAR_DKW(right), 
position_removed_from_DKW_set(right.position_removed_from_DKW_set),
position_removed_from_DKW(right.position_removed_from_DKW)
{}

bool CAR_DKW_nominal::SetParametersNotConsidered()
{
	if (dataP)
	{
		int offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols;
                position_removed_from_DKW = TIndex(offset, offset+dataP->bcf_vec.cols); // delta_bcf, delta_bcfLT
                offset = CAR_DKW::NumberParameters()-1;
                position_removed_from_DKW += TIndex(offset-dataP->ytips_vec.cols+1, offset);
                position_removed_from_DKW_set = true;
	}
	else 
		position_removed_from_DKW_set = false; 
	return position_removed_from_DKW_set; 
}

bool CAR_DKW_nominal::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (!position_removed_from_DKW_set && !SetParametersNotConsidered() )
		return false; 
	int DKW_n = CAR_DKW::NumberParameters(); 
	TDenseVector DKW_x(DKW_n, MINUS_INFINITY); 
	TIndex DKW_index(0,DKW_n-1); 
	DKW_index.UniqSubtract(position_removed_from_DKW); 
	DKW_x.Insert(DKW_index, _parameter, 0, DKW_n-1-position_removed_from_DKW.size); 
	return CAR_DKW::SetParameters_FromVectorToMatrix(DKW_x); 
}

bool CAR_DKW_nominal::SetParameters_InitializeABOmega()
{
	KAPPA_rn = KAPPA + SIGMAlambda1;
	try {
                KAPPA_inverse = Inverse(KAPPA);
                SIGMA_inverse = Inverse(SIGMA);
		Inv_KAPPA_rn = Inverse(KAPPA_rn);
                Inv_Kron_KAPPA_rn = Inverse(Kron(-1.0*KAPPA_rn, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA_rn));
        }
        catch (...) {
                return false;
        }

	TDenseVector KAPPAtheta = KAPPA * theta;
	if(!YieldFacLoad(ay, by, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, KAPPAtheta, rho0, rho1, lambda0, dataP->MATgrid))
                return false;
	TDenseMatrix Bx = MatrixExp(-dataP->Dt*KAPPA);
        TDenseVector Ax = (Identity(Nfac)-Bx)*theta;

        B.Zeros(Nfac+1, Nfac+1);
        B(0,0) = 1;
        B.InsertRowMatrix(0, 1, dataP->Dt*rho1_pi);
        B.Insert(1, 1, Bx);

        A.Zeros(Nfac+1);
        A(0) = rho0_pi*dataP->Dt;
        A.Insert(1, Ax);

        TDenseMatrix tmp_kron;
        try {
                tmp_kron = Inverse(Kron(-1.0*KAPPA, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA));
        }
        catch (...) {
                return false;
        }
        TDenseMatrix tmp_Omega = MultiplyTranspose(SIGMA, SIGMA);
        TDenseVector tmpVec(tmp_Omega.rows*tmp_Omega.cols);
        for (int i=0; i<tmp_Omega.cols; i++)
                tmpVec.Insert(i*tmp_Omega.rows, tmp_Omega.ColumnVector(i));
        tmpVec = -1.0*tmp_kron*tmpVec;
        TDenseMatrix OMEGA_x_ss(Nfac,Nfac);
        for (int i=0; i<OMEGA_x_ss.cols; i++)
                OMEGA_x_ss.InsertColumnMatrix(0, i, tmpVec.SubVector(i*OMEGA_x_ss.rows, (i+1)*OMEGA_x_ss.rows-1));

        tmp_Omega = Bx * SIGMA * Transpose(SIGMA) * Transpose(Bx) - SIGMA * Transpose(SIGMA);
        tmpVec.Zeros(tmp_Omega.rows*tmp_Omega.cols);
        for (int i=0; i<tmp_Omega.cols; i++)
                tmpVec.Insert(i*tmp_Omega.rows, tmp_Omega.ColumnVector(i));
        tmpVec = tmp_kron * tmpVec;
        TDenseMatrix OMEGA_x(Nfac,Nfac);
        for (int i=0; i<OMEGA_x.cols; i++)
                OMEGA_x.InsertColumnMatrix(0, i, tmpVec.SubVector(i*OMEGA_x.rows, (i+1)*OMEGA_x.rows-1));

        TDenseMatrix OMEGA_ss(Nfac+1, Nfac+1, 0.0);
        OMEGA_ss(0,0) = 1.0E6;
        OMEGA_ss.Insert(1,1,OMEGA_x_ss);

        double OMEGA_q=dataP->Dt*(InnerProduct(sigq,sigq) + sigqx*sigqx);
        TDenseVector OMEGA_xq;
        OMEGA_xq =Multiply(KAPPA_inverse, (Identity(Nfac)-Bx)*SIGMA*sigq);

        OMEGA.Zeros(Nfac+1,Nfac+1);
        OMEGA(0,0) = OMEGA_q;
        OMEGA.InsertRowMatrix(0,1,OMEGA_xq);
        OMEGA.InsertColumnMatrix(1,0,OMEGA_xq);
        OMEGA.Insert(1,1,OMEGA_x);

	int i=0;
        while (i<dataP->p_vec.rows && dataP->p_vec(i,0) == MINUS_INFINITY)
                i++;
        x0_0.Zeros(theta.Dimension()+1);
        x0_0(0) = dataP->p_vec(i,0);
        x0_0.Insert(1,theta);

        P0_0 = OMEGA_ss;
        return true;
}

TDenseVector CAR_DKW_nominal::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	if (dataP->p_vec(j,0) == MINUS_INFINITY)
        {
                b.Zeros(dataP->MATgrid.Dimension(), by.cols+1);
                b.Insert(0,1,by);
                a.CopyContent(ay);
                y_obs.CopyContent(dataP->y_vec.RowVector(j));
                R.Zeros(delta_y.Dimension(),delta_y.Dimension());
                for (int i=0; i<delta_y.Dimension(); i++)
                        R(i,i) = delta_y(i)*delta_y(i);
        }
        else
        {
                b.Zeros(dataP->MATgrid.Dimension()+1, Nfac+1);
                b(0,0) = 1;
                b.Insert(1,1,by);
                a.Zeros(ay.Dimension()+1);
                a.Insert(1,ay);
                y_obs.Zeros(dataP->p_vec.cols + dataP->y_vec.cols);
                y_obs.Insert(0, dataP->p_vec.RowVector(j));
                y_obs.Insert(dataP->p_vec.cols, dataP->y_vec.RowVector(j));
                R.Zeros(1+delta_y.Dimension(), 1+delta_y.Dimension());
                R(0,0) = delta_p * delta_p;
                for (int i=0; i<delta_y.Dimension(); i++)
                        R(i+1, i+1) = delta_y(i)*delta_y(i);
        }
	return TDenseVector(); 	
}

int CAR_DKW_nominal::NumberParameters() const
{
        int N = CAR_DKW::NumberParameters();
        N -= position_removed_from_DKW.size; // minus parameters of DKW that are not used 
        return N;
}

TDenseVector CAR_DKW_nominal::GetParameters()
{
	if (!if_internal_parameter_set)
	{
		TDenseVector complete_parameter(NumberParameters() + position_removed_from_DKW.size, 0.0); 
		complete_parameter.Insert(0, CAR_DKW::GetParameters()); 

		TIndex index_cleaned(0, NumberParameters()+position_removed_from_DKW.size-1); 
		index_cleaned.UniqSubtract(position_removed_from_DKW); 
		TDenseVector parameter_cleaned = complete_parameter.SubVector(index_cleaned); 
		internal_parameter.CopyContent(parameter_cleaned); 
		if_internal_parameter_set = true; 
	}
	return internal_parameter; 
}

bool CAR_DKW_nominal::SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
        if (iequals(which_one, std::string("delta_bcf")) )     
		return false; 
	else if (iequals(which_one, std::string("delta_tips")) )
		return false; 
	// because rho1pi and sigq are after delta_bcf and delta_bcfLT
	else if (iequals(which_one, std::string("rho1_pi")) )
	{
		int offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + 1; 
		TIndex local_fixed_index = FixedVariableParameter(rho1_pi, v, offset);	
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("sigq")) )
	{
		int offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + 1 + Nfac; 
		TIndex local_fixed_index = FixedVariableParameter(sigq, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else 
		return CAR_DKW::SetAsFixed(v, which_one); 
}

bool CAR_DKW_nominal::SetAsFixed(double v, const std::string &which_one)
{
	if (iequals(which_one, std::string("delta_bcfLT")) )
		return false; 
	// because rho0_pi, siqqx and delta_p are after delta_bcf and delta_bcfLT
	else if (iequals(which_one, std::string("rho0_pi")) )
	{
		rho0_pi = v; 
		TIndex local_fixed_index(Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols); 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("sigqx")) )
	{
		sigqx = v; 
		TIndex local_fixed_index(NumberParameters()-2); 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("delta_p")) )
	{
		delta_p = v; 
		TIndex local_fixed_index(NumberParameters()-1); 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else 
		return CAR_DKW::SetAsFixed(v, which_one); 
}
