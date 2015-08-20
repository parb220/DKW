#include "CAR_DKWl.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"
#include "prcsn.h"

CAR_DKWl::CAR_DKWl(int _Nfac, CData_FRBA *_dataP) : 
CAR_DKW(_Nfac, _dataP),
KAPPA_L(MINUS_INFINITY),
SIGMA_L(MINUS_INFINITY),
theta_L(MINUS_INFINITY),
rho1_L(_Nfac,MINUS_INFINITY),
rhoL_L(MINUS_INFINITY),
lambda0_L(MINUS_INFINITY),
SIGMAlambda1_L(MINUS_INFINITY), 
ay_TR(), by_TR(),
ay_L(), by_L()
{}

CAR_DKWl::CAR_DKWl(const CAR_DKWl &right) : 
CAR_DKW(right), 
KAPPA_L(right.KAPPA_L),
SIGMA_L(right.SIGMA_L),
theta_L(right.theta_L),
rho1_L(right.rho1_L),
rhoL_L(right.rhoL_L),
lambda0_L(right.lambda0_L),
SIGMAlambda1_L(right.SIGMAlambda1_L),
ay_TR(right.ay_TR),
by_TR(right.by_TR),
ay_L(right.ay_L),
by_L(right.by_L)
{
	rho1_L.CopyContent(right.rho1_L); 
	ay_TR.CopyContent(right.ay_TR); 
	by_TR.CopyContent(right.by_TR); 
	ay_L.CopyContent(right.ay_L); 
	by_L.CopyContent(right.by_L); 
}

TDenseVector CAR_DKWl::GetParameters()
{
	if (!if_internal_parameter_set)
	{
		TDenseVector complete_parameter(NumberParameters(),0.0);
                complete_parameter.Insert(0, CAR_DKW::GetParameters());
		int counter = CAR_DKW::NumberParameters(); 
		
		complete_parameter(counter) = KAPPA_L;
                counter ++;

		complete_parameter(counter) = SIGMA_L; 
                counter ++;

		complete_parameter(counter) = theta_L; 
		counter ++; 

		if (rho1_L.Dimension() != Nfac)
			rho1_L = TDenseVector(Nfac, MINUS_INFINITY); 
		complete_parameter.Insert(counter, rho1_L); 
		counter += Nfac;

		complete_parameter(counter) = rhoL_L; 
		counter ++; 

		complete_parameter(counter) = lambda0_L; 
		counter ++; 

		complete_parameter(counter) = SIGMAlambda1_L; 
		counter ++; 

		internal_parameter.CopyContent(complete_parameter);
                if_internal_parameter_set = true;
	}
	return internal_parameter; 
}

bool CAR_DKWl::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (_parameter.Dimension() < CAR_DKWl::NumberParameters())
		return false; 

	if (!CAR_DKW::SetParameters_FromVectorToMatrix(_parameter.SubVector(0,CAR_DKW::NumberParameters()-1)) )
		return false; 

	int counter = CAR_DKW::NumberParameters(); 
	KAPPA_L = _parameter(counter);
        counter ++;

	SIGMA_L = _parameter(counter); 
	counter ++; 

	theta_L = _parameter(counter); 
	counter ++; 

	if (rho1_L.Dimension() != Nfac)
		rho1_L.Resize(Nfac); 
	rho1_L.Insert(0,_parameter, counter, counter+Nfac-1); 
	counter += Nfac; 

	rhoL_L = _parameter(counter); 
	counter ++; 

	lambda0_L = _parameter(counter); 
	counter ++; 
	
	SIGMAlambda1_L = _parameter(counter); 
	counter ++; 
	return true; 
}

bool CAR_DKWl::SetParameters_InitializeABOmega()
{
	KAPPA_rn = KAPPA + SIGMAlambda1;
	try {
		KAPPA_inverse = Inverse(KAPPA); 
		SIGMA_inverse = Inverse(SIGMA); 
                Inv_KAPPA_rn = Inverse(KAPPA_rn);
                Inv_Kron_KAPPA_rn = Inverse(Kron(-1.0*KAPPA_rn, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA_rn));
	}
	catch(...) {
		return false; 
	}
	
	// ay, by
        TDenseVector KAPPAtheta = KAPPA * theta;
	if(!YieldFacLoad(ay, by, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, KAPPAtheta, rho0, rho1,lambda0, dataP->MATgrid))
                return false;

	// ay_R, by_R
	double rho0_R = rho0 - rho0_pi - 0.5*(InnerProduct(sigq, sigq) + sigqx*sigqx) + InnerProduct(lambda0, sigq);
	TDenseMatrix lambda1 = Multiply(SIGMA_inverse, SIGMAlambda1); 
        TDenseVector rho1_R = rho1 - rho1_pi + TransposeMultiply(lambda1, sigq);
        TDenseVector lambda0_R = lambda0 - sigq;
        if (!YieldFacLoad(ay_R, by_R, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, KAPPAtheta, rho0_R, rho1_R, lambda0_R, dataP->TIPSgrid))
                return false;

	// af, bf
	TDenseVector ForecastHor(2);
        ForecastHor[0] = 0.5; ForecastHor[1] = 1.0;
        af.Zeros(ForecastHor.Dimension());
        bf.Zeros(ForecastHor.Dimension(), Nfac);
        for (int i=0; i<ForecastHor.Dimension(); i++)
        {
                TDenseMatrix tempM = MatrixExp(-ForecastHor[i]*KAPPA);
                af(i) = ay(0)+InnerProduct(by.RowVector(0), theta, Identity(Nfac)-tempM);
                bf.InsertRowMatrix(i, 0, by.RowVector(0)*tempM);
        }

	// A B
	TDenseMatrix Bx = MatrixExp(-dataP->Dt*KAPPA);
        TDenseVector Ax = (Identity(Nfac)-Bx)*theta;

        double BL = exp(-dataP->Dt*KAPPA_L);
        double AL = (1.0-BL) * theta_L;

	B.Zeros(Nfac+2, Nfac+2);
        B(0,0) = 1.0;
        B.InsertRowMatrix(0,1,dataP->Dt*rho1_pi);
        B.Insert(1,1, Bx);
        B(Nfac+1, Nfac+1) = BL;

	A.Zeros(Nfac+2);
        A(0) = dataP->Dt*rho0_pi;
        A.Insert(1,Ax);
        A(Nfac+1) = AL;

	// OMEGA
	TDenseMatrix tmp_kron;
        try {
                tmp_kron = Inverse(Kron(-1.0*KAPPA, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA));
        }
        catch(...) {
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

	double OMEGA_L_ss = SIGMA_L*SIGMA_L/(2.0*KAPPA_L);
        double OMEGA_L=1.0/(2.0*KAPPA_L)*(1.0-exp(-2.0*KAPPA_L*dataP->Dt))*SIGMA_L*SIGMA_L;

	TDenseMatrix OMEGA_ss(Nfac+2, Nfac+2, 0.0);
        OMEGA_ss(0,0)=1.0e6;
        OMEGA_ss.Insert(1,1,OMEGA_x_ss);
        OMEGA_ss(Nfac+1,Nfac+1) = OMEGA_L_ss;

        double OMEGA_q = dataP->Dt*(InnerProduct(sigq, sigq) + sigqx*sigqx);
        TDenseVector OMEGA_xq;
        OMEGA_xq = Multiply(KAPPA_inverse, (Identity(Nfac)-Bx)*SIGMA*sigq);
        
	OMEGA.Zeros(Nfac+2,Nfac+2);
	OMEGA(0,0) = OMEGA_q;
        OMEGA.InsertRowMatrix(0,1,OMEGA_xq);
        OMEGA.InsertColumnMatrix(1,0,OMEGA_xq);
        OMEGA.Insert(1,1,OMEGA_x);
        OMEGA(Nfac+1,Nfac+1) = OMEGA_L;
	
	// x0_0, P0_0
	int i=0;
        while (i<dataP->p_vec.rows && dataP->p_vec(i,0) == MINUS_INFINITY)
                i++;
        x0_0.Zeros(theta.Dimension()+2);
        x0_0(0) = dataP->p_vec(i,0);
        x0_0.Insert(1,theta);
        x0_0(theta.Dimension()+1) = theta_L;

        P0_0 = OMEGA_ss; 

	// ay_TR, by_TR
	if (!YieldFacLoad(ay_TR, by_TR, KAPPA_rn, Inv_KAPPA_rn, Inv_Kron_KAPPA_rn, SIGMA, KAPPAtheta, rho0_R, rho1_R+rho1_L, lambda0_R,dataP->TIPSgrid))
                return false;

	// ay_L, by_L
	double KAPPA_L_rn = KAPPA_L + SIGMAlambda1_L;
	double Inv_KAPPA_L_rn = 1.0/(KAPPA_L_rn); 
	double Inv_Kron_KAPPA_L_rn = 1.0/((-1.0*KAPPA_L_rn)+(-1.0*KAPPA_L_rn)); 
	double KAPPA_Ltheta = KAPPA_L*theta_L; 

	if (!YieldFacLoad(ay_L, by_L, TDenseMatrix(1,1,KAPPA_L_rn), TDenseMatrix(1,1,Inv_KAPPA_L_rn), TDenseMatrix(1,1,Inv_Kron_KAPPA_L_rn), TDenseMatrix(1,1,SIGMA_L),TDenseVector(1,KAPPA_Ltheta), 0, TDenseVector(1,rhoL_L), TDenseVector(1,lambda0_L), dataP->TIPSgrid))
                return false;

        return true;
}

TDenseVector CAR_DKWl::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
        if (dataP->p_vec(j,0) == MINUS_INFINITY)
        {
                b.Zeros(dataP->MATgrid.Dimension(), by.cols+2);
                b.Insert(0, 1, by);
                a.CopyContent(ay); 
                y_obs.CopyContent(dataP->y_vec.RowVector(j)); 
                R.Zeros(delta_y.Dimension(), delta_y.Dimension());
                for (int i=0; i<delta_y.Dimension(); i++)
                        R(i,i) = delta_y(i) * delta_y(i);
        }
        else
        {
                b.Zeros(dataP->MATgrid.Dimension()+1, by.cols+2);
                b(0,0) = 1.0; 
                b.Insert(1,1,by); 
                a.Zeros(ay.Dimension()+1);
                a.Insert(1,ay); 
                y_obs.Zeros(dataP->p_vec.cols + dataP->y_vec.cols);
                y_obs.Insert(0, dataP->p_vec.RowVector(j)); 
                y_obs.Insert(dataP->p_vec.cols, dataP->y_vec.RowVector(j)); 
                R.Zeros(dataP->p_vec.cols+dataP->y_vec.cols, dataP->p_vec.cols+dataP->y_vec.cols);
                for (int i=0; i<dataP->p_vec.cols; i++)
                        R(i,i) = delta_p*delta_p;  
                for (int i=0; i<dataP->y_vec.cols; i++)
                        R(i+dataP->p_vec.cols, i+dataP->p_vec.cols) = delta_y(i) * delta_y(i);
        }

	if (dataP->bcf_vec(j,0) != MINUS_INFINITY)
        {
		//  b = [b;[zeros(size(bf,1),1),bf,zeros(size(bf,1),1)]];
		TDenseMatrix tmp_b(b);
                tmp_b.CopyContent(b);
                b.Zeros(tmp_b.rows+bf.rows, tmp_b.cols);
                b.Insert(0,0,tmp_b);
                b.Insert(tmp_b.rows, 1, bf);
                // a = [a;af];
                TDenseVector tmp_a(a);
                tmp_a.CopyContent(a);
                a.Zeros(tmp_a.Dimension()+af.Dimension());
                a.Insert(0,tmp_a);
                a.Insert(tmp_a.Dimension(), af);
                // y_obs= [y_obs; bcf_vec(j,:)'];
                TDenseVector tmp_y_obs(y_obs);
                tmp_y_obs.CopyContent(y_obs);
                y_obs.Zeros(y_obs.Dimension()+dataP->bcf_vec.cols);
                y_obs.Insert(0,tmp_y_obs);
                y_obs.Insert(tmp_y_obs.Dimension(), dataP->bcf_vec.RowVector(j));
		// Rvec= [Rvec; delta_bcf.^2];
		TDenseMatrix tmp_R(R);
                tmp_R.CopyContent(R);
                R.Zeros(tmp_R.rows+delta_bcf.Dimension(), tmp_R.cols+delta_bcf.Dimension());
                R.Insert(0,0,tmp_R);
                for (int i=0; i<delta_bcf.Dimension(); i++)
                        R(i+tmp_R.rows, i+tmp_R.cols) = delta_bcf(i)*delta_bcf(i);
	}

	if (dataP->bcfLT_vec(j,0) != MINUS_INFINITY)
        {
                double hor = dataP->horLT(j,0);
                TDenseMatrix W = 0.2* Multiply(KAPPA_inverse, (MatrixExp(-hor*KAPPA)-MatrixExp(-(hor+5.0)*KAPPA)) );
                double afLT = ay(0) + InnerProduct(by.RowVector(0), theta, Identity(Nfac)-W);
                TDenseVector bfLT = by.RowVector(0) * W;

                TDenseMatrix tmp_b(b);
                tmp_b.CopyContent(b);
                b.Zeros(tmp_b.rows+1, tmp_b.cols);
                b.Insert(0,0,tmp_b);
                b.InsertRowMatrix(tmp_b.rows, 1, bfLT);

                TDenseVector tmp_a(a);
                tmp_a.CopyContent(a);
                a.Zeros(tmp_a.Dimension()+1);
                a.Insert(0,tmp_a);
                a(tmp_a.Dimension()) = afLT;

                TDenseVector tmp_y_obs(y_obs);
                tmp_y_obs.CopyContent(y_obs);
                y_obs.Zeros(tmp_y_obs.Dimension()+1);
                y_obs.Insert(0,tmp_y_obs);
                y_obs(tmp_y_obs.Dimension()) = dataP->bcfLT_vec(j,0);

                TDenseMatrix tmp_R(R);
                tmp_R.CopyContent(R);
                R.Zeros(tmp_R.rows+1, tmp_R.cols+1);
                R.Insert(0,0,tmp_R);
                R(tmp_R.rows,tmp_R.cols) = delta_bcfLT * delta_bcfLT;
        }
	
	if (dataP->ytips_vec(j,0) != MINUS_INFINITY)
        {
                TDenseMatrix tmp_b(b);
                tmp_b.CopyContent(b);
                b.Zeros(tmp_b.rows+by_TR.rows, tmp_b.cols);
                b.Insert(0,0,tmp_b);
                b.Insert(tmp_b.rows,1,by_TR);
                b.Insert(tmp_b.rows,1+by_TR.cols, by_L);

                TDenseVector tmp_a(a);
                tmp_a.CopyContent(a);
                a.Zeros(tmp_a.Dimension()+ay_L.Dimension());
                a.Insert(0, tmp_a);
                a.Insert(tmp_a.Dimension(), ay_L+ay_TR);

                TDenseVector tmp_y_obs(y_obs);
                tmp_y_obs.CopyContent(y_obs);
                y_obs.Zeros(tmp_y_obs.Dimension()+dataP->ytips_vec.cols);
                y_obs.Insert(0,tmp_y_obs);
                y_obs.Insert(tmp_y_obs.Dimension(), dataP->ytips_vec.RowVector(j));

                TDenseMatrix tmp_R(R);
                tmp_R.CopyContent(R);
                R.Zeros(tmp_R.rows+delta_tips.Dimension(), tmp_R.cols+delta_tips.Dimension());
                R.Insert(0,0,tmp_R);
                for (int i=0; i<delta_tips.Dimension(); i++)
                        R(i+tmp_R.rows, i+tmp_R.cols) = delta_tips(i)*delta_tips(i);
        }
	return TDenseVector(); 
}

bool CAR_DKWl::SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
        TIndex local_fixed_index;
        int offset;
        if (iequals(which_one, std::string("rho1_L")) )
        {
                offset = CAR_DKW::NumberParameters() +  1 + 1 + 1;
                local_fixed_index = FixedVariableParameter(rho1_L, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else
                return CAR_DKW::SetAsFixed(v, which_one);
}

bool CAR_DKWl::SetAsFixed(double v, const std::string &which_one)
{
	TIndex local_fixed_index;
        if (iequals(which_one, std::string("KAPPA_L"))  )
        {
                KAPPA_L = v;
                local_fixed_index += CAR_DKW::NumberParameters(); 
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if (iequals(which_one, std::string("SIGMA_L"))  )
        {
                SIGMA_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() +1; 
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if (iequals(which_one, std::string("theta_L")) )
        {
                theta_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() + 1 + 1;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if(iequals(which_one, std::string("rhoL_L")) )
        {
                rhoL_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() + 1 + 1 + 1 + Nfac;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if(iequals(which_one, std::string("lambda0_L")) )
        {
                lambda0_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() + 1 + 1 + 1 + Nfac + 1;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if(iequals(which_one, std::string("SIGMAlambda1_L")) )
        {
                SIGMAlambda1_L = v;  
                local_fixed_index += CAR_DKW::NumberParameters() + 1 + 1 + 1 + Nfac + 1 + 1;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else
                return CAR_DKW::SetAsFixed(v, which_one);
}

int CAR_DKWl::NumberParameters() const
{
        return CAR_DKW::NumberParameters() + Nfac + 6;
}


