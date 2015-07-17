#include "prcsn.h"
#include "CAR_DKW.hpp"
#include "CData_FRBA.hpp"

TDenseMatrix MatrixExp(const TDenseMatrix &M); 

CAR_DKW::CAR_DKW(int _Nfac, CData_FRBA* _dataP) : 
CAR(_Nfac, _dataP),
KAPPA(TDenseMatrix(_Nfac,_Nfac,MINUS_INFINITY)), SIGMA(TDenseMatrix(_Nfac,_Nfac,MINUS_INFINITY)), theta(TDenseVector(_Nfac,MINUS_INFINITY)), rho0(MINUS_INFINITY), rho1(TDenseVector(_Nfac,MINUS_INFINITY)), 
lambda0(TDenseVector(_Nfac,MINUS_INFINITY)), SIGMAlambda1(TDenseMatrix(_Nfac,_Nfac,MINUS_INFINITY)), delta_y(), delta_bcf(),
delta_bcfLT(MINUS_INFINITY), rho0_pi(MINUS_INFINITY), rho1_pi(TDenseVector(_Nfac,MINUS_INFINITY)), sigq(TDenseVector(_Nfac,MINUS_INFINITY)), sigqx(MINUS_INFINITY),
delta_p(MINUS_INFINITY), delta_tips(), 
ay(), by(), ay_R(), by_R(),af(),bf(),
FixedIndex(), VariableIndex(), if_variable_index_set(false),
internal_parameter(), if_internal_parameter_set(false)
{
	if(dataP)
	{
		delta_y = TDenseVector(dataP->y_vec.cols, MINUS_INFINITY); 
		delta_bcf = TDenseVector(dataP->bcf_vec.cols, MINUS_INFINITY); 
		delta_tips = TDenseVector(dataP->ytips_vec.cols, MINUS_INFINITY); 
	}
}

CAR_DKW::CAR_DKW(const CAR_DKW &right) :
CAR(right),
KAPPA(right.KAPPA), SIGMA(right.SIGMA), theta(right.theta), 
rho0(right.rho0), rho1(right.rho1), 
lambda0(right.lambda0), SIGMAlambda1(right.SIGMAlambda1),
delta_y(right.delta_y), delta_bcf(right.delta_bcf), delta_bcfLT(right.delta_bcfLT),
rho0_pi(right.rho0_pi), rho1_pi(right.rho1_pi), 
sigq(right.sigq), sigqx(right.sigqx), delta_p(right.delta_p), delta_tips(right.delta_tips), 
ay(right.ay), by(right.by), ay_R(right.ay_R), by_R(right.by_R),af(right.af),bf(right.bf),
FixedIndex(right.FixedIndex), VariableIndex(right.VariableIndex), if_variable_index_set(right.if_variable_index_set),
internal_parameter(right.internal_parameter), if_internal_parameter_set(right.if_internal_parameter_set)
{
	KAPPA.CopyContent(right.KAPPA); 
	SIGMA.CopyContent(right.SIGMA); 
	theta.CopyContent(right.theta); 
	rho1.CopyContent(right.rho1); 
	lambda0.CopyContent(right.lambda0); 
	SIGMAlambda1.CopyContent(right.SIGMAlambda1); 
	delta_y.CopyContent(right.delta_y); 
	delta_bcf.CopyContent(right.delta_bcf); 
	rho1_pi.CopyContent(right.rho1_pi); 
	sigq.CopyContent(right.sigq); 
	delta_tips.CopyContent(right.delta_tips); 
	ay.CopyContent(right.ay); 
	by.CopyContent(right.by); 
	ay_R.CopyContent(right.ay_R); 
	by_R.CopyContent(right.by_R); 
	af.CopyContent(right.af); 
	bf.CopyContent(right.bf); 
	internal_parameter.CopyContent(right.internal_parameter);
}

TDenseVector CAR_DKW::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	// p_vec
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

	// bcf
	if (dataP->bcf_vec(j,0) != MINUS_INFINITY)
	{
		// b = [b; [zeros(size(bf,1),1), bf]]
		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+bf.rows, tmp_b.cols); 
		b.Insert(0, 0, tmp_b); 
		b.Insert(tmp_b.rows, 1, bf); 
		// a = [a; af]
		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+af.Dimension()); 
		a.Insert(0, tmp_a); 
		a.Insert(tmp_a.Dimension(), af); 
		// y_obs = [[y_obs; bcf_vec(j,:)'];
		TDenseVector tmp_y_obs(y_obs); 
		tmp_y_obs.CopyContent(y_obs); 
		y_obs.Zeros(tmp_y_obs.Dimension() + dataP->bcf_vec.cols); 
		y_obs.Insert(0, tmp_y_obs); 
		y_obs.Insert(tmp_y_obs.Dimension(), dataP->bcf_vec.RowVector(j)); 
		// Rvec= [Rvec; delta_bcf.^2];
		TDenseMatrix tmp_R(R); 
		tmp_R.CopyContent(R); 
		R.Zeros(tmp_R.rows+delta_bcf.Dimension(), tmp_R.rows+delta_bcf.Dimension()); 
		R.Insert(0, 0, tmp_R); 
		for (int i=tmp_R.rows; i<tmp_R.rows+delta_bcf.Dimension(); i++)
			R(i,i) = delta_bcf(i-tmp_R.rows)*delta_bcf(i-tmp_R.rows); 
	}

	// bcfLT
	if (dataP->bcfLT_vec(j,0) != MINUS_INFINITY)
	{
		TDenseMatrix W = Inverse(KAPPA) * (MatrixExp(-dataP->horLT(j,0)*KAPPA)-MatrixExp(-(dataP->horLT(j,0)+5.0)*KAPPA)) * (1.0/5.0); 
		double afLT = ay(0) + InnerProduct(by.RowVector(0), theta, Identity(Nfac)-W ); 
		TDenseVector bfLT = by.RowVector(0) * W; 
		
		// b = [b;[0,bfLT]]
		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+1, tmp_b.cols); 
		b.Insert(0, 0, tmp_b); 
		b.InsertRowMatrix(tmp_b.rows, 1, bfLT); 
		// a = [a; afLT]
		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(a.Dimension()+1); 
		a.Insert(0, tmp_a); 
		a(tmp_a.Dimension()) = afLT; 
		// y_obs= [y_obs; bcfLT_vec(j)]
		TDenseVector tmp_y_obs(y_obs); 
		tmp_y_obs.CopyContent(y_obs); 
		y_obs.Zeros(tmp_y_obs.Dimension()+1); 
		y_obs.Insert(0,tmp_y_obs); 
		y_obs(tmp_y_obs.Dimension()) = dataP->bcfLT_vec(j,0); 
		// Rvec= [Rvec; delta_bcfLT^2];
		TDenseMatrix tmp_R(R); 
		tmp_R.CopyContent(R); 
		R.Zeros(tmp_R.rows+1, tmp_R.cols+1); 
		R.Insert(0, 0, tmp_R); 
		R(tmp_R.rows, tmp_R.cols) = delta_bcfLT *delta_bcfLT; 
	}

	// tips
	if (dataP->ytips_vec(j,0) != MINUS_INFINITY)
	{
		// b = [b;[zeros(size(by_R,1),1),by_R]];
		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+by_R.rows, tmp_b.cols); 
		b.Insert(0, 0, tmp_b); 
		b.Insert(tmp_b.rows, 1, by_R); 
				
		// a  = [a; ay_R]; 
		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+ay_R.Dimension()); 
		a.Insert(0,tmp_a); 
		a.Insert(tmp_a.Dimension(), ay_R); 

		// y_obs= [y_obs; ytips_vec(j,:)'];
		TDenseVector tmp_y_obs(y_obs); 
		tmp_y_obs.CopyContent(y_obs); 
		y_obs.Zeros(tmp_y_obs.Dimension()+dataP->ytips_vec.cols); 
		y_obs.Insert(0, tmp_y_obs); 
		y_obs.Insert(tmp_y_obs.Dimension(), dataP->ytips_vec.RowVector(j)); 
	
		// Rvec= [Rvec; delta_tips.^2]
		TDenseMatrix tmp_R(R); 
		tmp_R.CopyContent(R); 
		R.Zeros(tmp_R.rows+delta_tips.Dimension(), tmp_R.cols+delta_tips.Dimension()); 
		R.Insert(0, 0, tmp_R); 
		for (int i=tmp_R.rows; i<tmp_R.rows+delta_tips.Dimension(); i++)
			R(i,i) = delta_tips(i-tmp_R.rows)*delta_tips(i-tmp_R.rows); 
	}
	return TDenseVector(0); 
}

int CAR_DKW::NumberParameters() const
{
	int n=0; 
	n += Nfac*Nfac; // KAPPA
	n += Nfac*Nfac; // SIGMA
	n += Nfac; 	// theta
	n ++; 		// rho0
	n += Nfac;	// rho1
	n += Nfac;	// lambda0
	n += Nfac*Nfac;	// SIGMAlambda1
	n += dataP->y_vec.cols;	// delta_y
	n += dataP->bcf_vec.cols; // delta_bcf; 
	n += dataP->bcfLT_vec.cols;	// delta_bcfLT

	n ++; 		// rho0_pi
	n += Nfac; 	// rho1_pi;
	n += Nfac; 	// sigq; 
	n ++; 		// sigqx; 
	
	n ++;		// delta_p;
	n += dataP->ytips_vec.cols;	// delta_tip

	return n; 
}

int CAR_DKW::NumberVariableParameters() const
{
	return NumberParameters()-FixedIndex.size; 
}

int CAR_DKW::NumberFixedParameters() const
{
	return FixedIndex.size; 
}

TIndex CAR_DKW::GetVariableIndex()
{
	if (!if_variable_index_set)
	{
		TIndex all(0, NumberParameters()-1); 
		VariableIndex = all.UniqSubtract(FixedIndex); 
		if_variable_index_set = true; 
	}	
	return VariableIndex; 
}

bool CAR_DKW::SetParameters(const TDenseVector &_parameter)
{
	TDenseVector current_parameter = GetParameters(); 
	if (_parameter.Dimension() == NumberVariableParameters())
		current_parameter.Insert(GetVariableIndex(), _parameter); 	
	else if (_parameter.Dimension() == NumberParameters())
		current_parameter.Insert(GetVariableIndex(), _parameter(GetVariableIndex()));
	else 
		return false; 
	if (!SetParameters_FromVectorToMatrix(current_parameter))
		return false; 
	if (!SetParameters_ABOmega())
		return false; 
	internal_parameter = current_parameter;
	if_internal_parameter_set = true; 
	return true; 
}

bool CAR_DKW::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (_parameter.Dimension() < CAR_DKW::NumberParameters())
		return false; 

	int counter = 0; 
	if (KAPPA.rows != Nfac || KAPPA.cols != Nfac)
		KAPPA.Resize(Nfac, Nfac);
	for (int i=0; i<KAPPA.cols; i++)
	{
		KAPPA.InsertColumnMatrix(0, i, _parameter, counter, counter+KAPPA.rows-1); 
		counter += KAPPA.rows; 
	}

	if (SIGMA.rows != Nfac || SIGMA.cols != Nfac)
		SIGMA.Resize(Nfac, Nfac); 
	for (int i=0; i<SIGMA.cols; i++)
	{
		SIGMA.InsertColumnMatrix(0, i, _parameter, counter, counter+SIGMA.rows-1); 
		counter += SIGMA.rows; 
	}

	if (theta.Dimension() != Nfac)
		theta.Resize(Nfac); 
	theta.Insert(0, _parameter, counter, counter+theta.Dimension()-1); 
	counter += theta.Dimension(); 

	rho0 = _parameter(counter); 
	counter ++; 

	if (rho1.Dimension() != Nfac)
		rho1.Resize(Nfac); 
	rho1.Insert(0, _parameter, counter, counter+rho1.Dimension()-1); 
	counter += rho1.Dimension(); 

	if (lambda0.Dimension() != Nfac)
		lambda0.Resize(Nfac); 
	lambda0.Insert(0, _parameter, counter, counter+lambda0.Dimension()-1); 
	counter += lambda0.Dimension(); 

	if (SIGMAlambda1.rows != Nfac || SIGMAlambda1.cols != Nfac)
		SIGMAlambda1.Resize(Nfac, Nfac); 
	for (int i=0; i<SIGMAlambda1.cols; i++)
	{
		SIGMAlambda1.InsertColumnMatrix(0, i, _parameter, counter, counter+SIGMAlambda1.rows-1); 
		counter += SIGMAlambda1.rows; 
	}

	if (delta_y.Dimension() != dataP->y_vec.cols)
		delta_y.Resize(dataP->y_vec.cols); 
	delta_y.Insert(0, _parameter, counter, counter+delta_y.Dimension()-1); 
	counter += delta_y.Dimension();

	if (delta_bcf.Dimension() != dataP->bcf_vec.cols)
		delta_bcf.Resize(dataP->bcf_vec.cols); 
	delta_bcf.Insert(0, _parameter, counter, counter+delta_bcf.Dimension()-1); 
	counter += delta_bcf.Dimension(); 

	delta_bcfLT =  _parameter(counter); 
	counter ++; 

	rho0_pi = _parameter(counter); 
	counter ++; 

	if (rho1_pi.Dimension() != Nfac)
		rho1_pi.Resize(Nfac); 
	rho1_pi.Insert(0, _parameter, counter, counter+rho1_pi.Dimension()-1); 
	counter += rho1_pi.Dimension(); 

	if (sigq.Dimension() != Nfac)	
		sigq.Resize(Nfac); 
	sigq.Insert(0, _parameter, counter, counter+sigq.Dimension()-1);
	counter += sigq.Dimension(); 

	sigqx = _parameter(counter); 
	counter++; 

	delta_p = _parameter(counter); 
	counter ++; 

	if (delta_tips.Dimension() != dataP->ytips_vec.cols)
		delta_tips.Resize(dataP->ytips_vec.cols); 
	delta_tips.Insert(0, _parameter, counter, counter+delta_tips.Dimension()-1); 
	counter += delta_tips.Dimension(); 
	return true; 
}

bool CAR_DKW::SetParameters_ABOmega()
{
	TDenseMatrix lambda1 = InverseMultiply(SIGMA, SIGMAlambda1);
        YieldFacLoad(ay, by, KAPPA, SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,dataP->MATgrid);

        double rho0_R = rho0 - rho0_pi - 0.5*(InnerProduct(sigq, sigq) + sigqx*sigqx) + InnerProduct(lambda0, sigq);
        TDenseVector rho1_R = rho1 - rho1_pi + TransposeMultiply(lambda1, sigq);
        TDenseVector lambda0_R = lambda0 - sigq;
        TDenseMatrix SIGMAlambda1_R = SIGMAlambda1;
        YieldFacLoad(ay_R, by_R, KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,dataP->TIPSgrid);

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

	TDenseMatrix Bx = MatrixExp(-dataP->Dt*KAPPA); 
	TDenseVector Ax = (Identity(Nfac)-Bx)*theta; 

	B.Zeros(Nfac+1, Nfac+1); 
	B(0,0) = 1; 
	B.InsertRowMatrix(0, 1, dataP->Dt*rho1_pi); 
	B.Insert(1, 1, Bx);

	A.Zeros(Nfac+1); 
	A(0) = rho0_pi*dataP->Dt;
	A.Insert(1, Ax); 

	TDenseMatrix tmp_kron = Inverse(Kron(-1.0*KAPPA, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA)); 
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
	TDenseVector OMEGA_xq=Inverse(KAPPA)*(Identity(Nfac)-Bx)*SIGMA*sigq;
	OMEGA.Zeros(Nfac+1,Nfac+1); 
	OMEGA(0,0) = OMEGA_q; 
	OMEGA.InsertRowMatrix(0,1,OMEGA_xq); 
	OMEGA.InsertColumnMatrix(1,0,OMEGA_xq); 
	OMEGA.Insert(1,1,OMEGA_x);

	// x0_0, P0_0 
	int i=0; 
	while (i<dataP->p_vec.rows && dataP->p_vec(i,0) == MINUS_INFINITY)
		i++; 
	x0_0.Zeros(theta.Dimension()+1); 
	x0_0(0) = dataP->p_vec(i,0); 
	x0_0.Insert(1,theta); 
	
	P0_0 = OMEGA_ss; 	
	return true; 
}	

TDenseVector CAR_DKW::GetParameters() 
{
	if (!if_internal_parameter_set)
	{
		int counter =0; 
		internal_parameter.Zeros(NumberParameters()); 
		// KAPPA
		if (KAPPA.rows != Nfac || KAPPA.cols != Nfac)
		{
			KAPPA = TDenseMatrix(Nfac, Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, Ones(Nfac*Nfac)*MINUS_INFINITY); 
		}
		else 
		{
			for (int i=0; i<KAPPA.cols; i++)
				internal_parameter.Insert(counter+i*KAPPA.rows, KAPPA.ColumnVector(i)); 
		}
		counter += Nfac * Nfac; 

		// SIGMA
		if (SIGMA.rows != Nfac || SIGMA.cols != Nfac)
		{
			SIGMA = TDenseMatrix(Nfac, Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, Ones(Nfac*Nfac)*MINUS_INFINITY);
		}
		else 
		{
			for (int i=0; i<SIGMA.cols; i++)
				internal_parameter.Insert(counter+i*SIGMA.rows, SIGMA.ColumnVector(i)); 
		}
		counter += Nfac * Nfac; 	

		// theta
		if (theta.Dimension() != Nfac)
		{
			theta = TDenseVector(Nfac, MINUS_INFINITY);
			internal_parameter.Insert(counter, Ones(Nfac)*MINUS_INFINITY); 
		} 
		else 
			internal_parameter.Insert(counter, theta); 
		counter += Nfac; 

		// rho0
		internal_parameter(counter) = rho0; 
		counter ++; 

		// rho1
		if (rho1.Dimension() != Nfac)
		{
			rho1 = TDenseVector(Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, Ones(Nfac)*MINUS_INFINITY); 
		}
		else 
			internal_parameter.Insert(counter, rho1); 
		counter += Nfac; 

		// lambda0
		if (lambda0.Dimension() != Nfac)
		{
			lambda0 = TDenseVector(Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, Ones(Nfac)*MINUS_INFINITY);
		}
		else 
			internal_parameter.Insert(counter, lambda0); 
		counter += Nfac; 

		// SIGMAlambda1
		if (SIGMAlambda1.rows != Nfac || SIGMAlambda1.cols != Nfac)
		{
			SIGMAlambda1 = TDenseMatrix(Nfac,Nfac,MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones(Nfac*Nfac)); 
		}
		else 
		{
			for (int i=0; i<SIGMAlambda1.cols; i++)
				internal_parameter.Insert(counter+i*SIGMAlambda1.rows, SIGMAlambda1.ColumnVector(i)); 
		}
		counter += Nfac * Nfac; 

		// delta_y
		if (delta_y.Dimension() != dataP->y_vec.cols)
		{
			delta_y = TDenseVector(dataP->y_vec.cols, MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones(dataP->y_vec.cols)); 
		}
		else 
			internal_parameter.Insert(counter, delta_y); 
		counter += dataP->y_vec.cols; 

		// delta_bcf
		if (delta_bcf.Dimension()  != dataP->bcf_vec.cols)
		{
			delta_bcf = TDenseVector(dataP->bcf_vec.cols, MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones(dataP->bcf_vec.cols)); 
		}
		else 
			internal_parameter.Insert(counter, delta_bcf); 
		counter += dataP->bcf_vec.cols; 

		// delta_bcfLT
		internal_parameter(counter) = delta_bcfLT; 
		counter ++;
 
		// rho0_pi
		internal_parameter(counter) = rho0_pi; 
		counter ++;
 
		// rho1_pi
		if (rho1_pi.Dimension() != Nfac)
		{
			rho1_pi = TDenseVector(Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones(Nfac)); 
		}
		else
			internal_parameter.Insert(counter, rho1_pi); 
		counter += Nfac;
 
		// sigq; 
		if (sigq.Dimension() != Nfac)
		{
			sigq = TDenseVector(Nfac, MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones(Nfac)); 
		}
		else 
			internal_parameter.Insert(counter, sigq); 
		counter += Nfac; 

		// sigqx
		internal_parameter(counter) = sigqx; 
		counter ++ ; 

		// delta_p
		internal_parameter(counter) = delta_p; 
		counter ++; 

		// delta_tips
		if (delta_tips.Dimension() != dataP->ytips_vec.cols)
		{
			delta_tips = TDenseVector( dataP->ytips_vec.cols, MINUS_INFINITY); 
			internal_parameter.Insert(counter, MINUS_INFINITY*Ones( dataP->ytips_vec.cols) ); 
		}
		else 
			internal_parameter.Insert(counter, delta_tips); 
		counter +=  dataP->ytips_vec.cols; 

		if_internal_parameter_set = true; 
	}
	return internal_parameter; 
}

bool CAR_DKW::SetAsFixed(const TDenseMatrix &M, const std::string &which_one)
{
	TIndex local_fixed_index; 
	if (iequals(which_one, std::string("KAPPA")) )
		local_fixed_index = FixedVariableParameter(KAPPA, M, 0); 
	else if(iequals(which_one, std::string("SIGMA")) )
	{
		int offset = Nfac * Nfac; 
		local_fixed_index = FixedVariableParameter(SIGMA, M, offset); 
	}
	else if (iequals(which_one, std::string("SIGMAlambda1")) )
	{
		int offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac; 
		local_fixed_index = FixedVariableParameter(SIGMAlambda1, M, offset);  	
	}
	else 
		return false; 

	FixedIndex.UniqMerge(local_fixed_index); 
	if_variable_index_set = false; 

	return true; 
}

bool CAR_DKW::SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
	TIndex local_fixed_index; 
	int offset; 
	if (iequals(which_one, std::string("theta")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac; 
		local_fixed_index = FixedVariableParameter(theta, v, offset); 	
	}
	else if (iequals(which_one, std::string("rho1")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1; 
		local_fixed_index = FixedVariableParameter(rho1, v, offset); 
	}
	else if (iequals(which_one, std::string("lambda0")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac; 
		local_fixed_index = FixedVariableParameter(lambda0, v, offset); 
	}
	else if (iequals(which_one, std::string("delta_y")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac; 
		local_fixed_index = FixedVariableParameter(delta_y, v, offset); 
	}
	else if (iequals(which_one, std::string("delta_bcf")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols; 
		local_fixed_index = FixedVariableParameter(delta_bcf, v, offset); 
	}
	else if (iequals(which_one, std::string("rho1_pi")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1; 
		local_fixed_index = FixedVariableParameter(rho1_pi, v, offset); 
	}
	else if (iequals(which_one, std::string("sigq")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac; 
		local_fixed_index = FixedVariableParameter(sigq, v, offset);
	}
	else if (iequals(which_one, std::string("delta_tips")) )
	{
		offset = Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac + Nfac + 1 + 1; 
		local_fixed_index = FixedVariableParameter(delta_tips, v, offset);
	}
	else 
		return false; 

	FixedIndex.UniqMerge(local_fixed_index); 
	if_variable_index_set = false; 
	return true; 
}

bool CAR_DKW::SetAsFixed(double v, const std::string &which_one)
{
	TIndex local_fixed_index; 
	if (iequals(which_one, std::string("rho0"))  )
	{
		rho0 = v; 
		local_fixed_index += (Nfac*Nfac + Nfac*Nfac + Nfac); 
	}
	else if (iequals(which_one, std::string("delta_bcfLT")) )
	{
		delta_bcfLT = v; 
		local_fixed_index += (Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols); 
	}
	else if (iequals(which_one, std::string("rho0_pi")) )
	{
		rho0_pi = v; 
		local_fixed_index += (Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 ); 
	}
	else if (iequals(which_one, std::string("sigqx")) )
	{
		sigqx = v; 
		local_fixed_index += (Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac + Nfac); 
	}
	else if (iequals(which_one, std::string("delta_p")) )
	{
		delta_p = v; 
		local_fixed_index += (Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac + Nfac + 1); 
	}
	else 
		return false; 

	FixedIndex.UniqMerge(local_fixed_index);   
        if_variable_index_set = false;	
	return true; 
}
