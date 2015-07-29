#include "CAR_DKWlv_o.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"
#include "prcsn.h"

CAR_DKWlv_o::CAR_DKWlv_o(int _Nfac, CData_FRBA *_dataP) : 
CAR_DKW(_Nfac, _dataP),
position_removed_from_DKW_set(false), 
position_removed_from_DKW(),   
rhov(MINUS_INFINITY), rhov_pi(MINUS_INFINITY), KAPPAv(MINUS_INFINITY), SIGMAv(MINUS_INFINITY), thetav(MINUS_INFINITY), 
gamv(MINUS_INFINITY), gamvx(MINUS_INFINITY), rho(MINUS_INFINITY), 
delta_options(), 
KAPPA_L(MINUS_INFINITY), SIGMA_L(MINUS_INFINITY), theta_L(MINUS_INFINITY), 
rho1_L(_Nfac, MINUS_INFINITY), 
rhoL_L(MINUS_INFINITY), lambda0_L(MINUS_INFINITY), SIGMAlambda1_L(MINUS_INFINITY), 
cy(), cy_R(), ay_TR(), by_TR(), cy_TR(), ay_L(), by_L(), cf(), aI_Q(), bI_Q(), cI_Q()
{
	if (dataP)
	{
		delta_options = TDenseVector(dataP->MATgrid_options.Dimension(), MINUS_INFINITY); 
		position_removed_from_DKW = TIndex(Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac + Nfac); 
		position_removed_from_DKW_set = true; 
	}
}

CAR_DKWlv_o::CAR_DKWlv_o(const CAR_DKWlv_o &right) : 
CAR_DKW(right), 
position_removed_from_DKW_set(right.position_removed_from_DKW_set),
position_removed_from_DKW(right.position_removed_from_DKW), 
rhov(right.rhov), rhov_pi(right.rhov_pi),  KAPPAv(right.KAPPAv), SIGMAv(right.SIGMAv), thetav(right.thetav), 
gamv(right.gamv), gamvx(right.gamvx), rho(right.rho), 
delta_options(right.delta_options), 
KAPPA_L(right.KAPPA_L), SIGMA_L(right.SIGMA_L), theta_L(right.theta_L), 
rho1_L(right.rho1_L), 
rhoL_L(right.rhoL_L), lambda0_L(right.lambda0_L), SIGMAlambda1_L(right.SIGMAlambda1_L), 
cy(right.cy), cy_R(right.cy_R), ay_TR(right.ay_TR), by_TR(right.by_TR), cy_TR(right.cy_TR),
ay_L(right.ay_L), by_L(right.by_L), cf(right.cf), aI_Q(right.aI_Q), bI_Q(right.bI_Q), cI_Q(right.cI_Q) 
{
	delta_options.CopyContent(right.delta_options); 
	rho1_L.CopyContent(right.rho1_L); 
	cy.CopyContent(right.cy);
	cy_R.CopyContent(right.cy_R);
	ay_TR.CopyContent(right.ay_TR);
	by_TR.CopyContent(right.by_TR);
	cy_TR.CopyContent(right.cy_TR);
	ay_L.CopyContent(right.ay_L);
	by_L.CopyContent(right.by_L);
	cf.CopyContent(right.cf);
	aI_Q.CopyContent(right.aI_Q);
	bI_Q.CopyContent(right.bI_Q);
	cI_Q.CopyContent(right.cI_Q);
}

bool CAR_DKWlv_o::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (!position_removed_from_DKW_set && !SetParametersNotConsidered())
		return false; 
	// for parameters of DKW
	int DKW_n = CAR_DKW::NumberParameters(); 
	TDenseVector DKW_x(DKW_n); 
	TIndex DKW_index(0, DKW_n-1); 
	DKW_index.UniqSubtract(position_removed_from_DKW); 
	DKW_x.Insert(DKW_index, _parameter, 0, DKW_n-1-position_removed_from_DKW.size); 
	for (int i=0; i<position_removed_from_DKW.size; i++)
		DKW_x(position_removed_from_DKW[i]) = MINUS_INFINITY; 
	if (!CAR_DKW::SetParameters_FromVectorToMatrix(DKW_x)) 
		return false; 
	
	// DKWlv_o specific parameters
	int counter = DKW_n - position_removed_from_DKW.size; 
	rhov = _parameter(counter); 
	counter ++; 

	rhov_pi = _parameter(counter); 
	counter ++; 

	KAPPAv = _parameter(counter); 
	counter ++; 

	SIGMAv = _parameter(counter); 
	counter ++; 

	thetav = _parameter(counter); 
	counter ++; 

	gamv = _parameter(counter); 
	counter ++; 

	gamvx = _parameter(counter); 
	counter ++; 

	rho = _parameter(counter); 
	counter ++; 

	if (delta_options.Dimension() != dataP->MATgrid_options.Dimension())
		delta_options.Resize(dataP->MATgrid_options.Dimension()); 
	delta_options.Insert(0, _parameter, counter, counter+delta_options.Dimension()-1); 
	counter += delta_options.Dimension(); 

	KAPPA_L = _parameter(counter); 
	counter ++; 

	SIGMA_L = _parameter(counter); 
	counter ++; 

	theta_L = _parameter(counter); 
	counter ++; 

	if (rho1_L.Dimension() != Nfac)
		rho1_L.Resize(Nfac); 
	rho1_L.Insert(0, _parameter, counter, counter+rho1_L.Dimension()-1); 
	counter += rho1_L.Dimension(); 

	rhoL_L = _parameter(counter); 
	counter ++; 

	lambda0_L = _parameter(counter); 
	counter ++; 

	SIGMAlambda1_L = _parameter(counter); 
	counter ++; 
	return true; 
}

bool CAR_DKWlv_o::SetParameters_InitializeABOmega()
{
	TDenseMatrix lambda1;
        try {
                lambda1 = InverseMultiply(SIGMA, SIGMAlambda1);
        }
        catch(...){
                return false;
        }

	TDenseMatrix SIGMA_SIGMA =  MultiplyTranspose(SIGMA,SIGMA); 
	// ay, by, cy
	TDenseMatrix KAPPA_rn = KAPPA + SIGMAlambda1; 
	TDenseVector KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0; 
	double KAPPAv_rn = KAPPAv+sqrt(1.0-rho*rho)*SIGMAv*gamv+rho*SIGMAv*gamvx;
	if (!YieldFacLoad_ODE3(ay,by,cy, dataP->MATgrid,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA_SIGMA,rho0,rho1,rhov,KAPPAv_rn,SIGMAv*SIGMAv))
		return false; 
	for (int i=0; i<ay.Dimension(); i++)
		ay(i) = -ay(i) / dataP->MATgrid(i); 
	TDenseMatrix X = RepRowVector(dataP->MATgrid, Nfac, 1); 
	for (int j=0; j<by.cols; j++)
		for (int i=0; i<by.rows; i++)
			by(i,j) = -by(i,j) / X(i,j); 
	by = Transpose(by); 
	for (int i=0; i<cy.Dimension(); i++)
		cy(i) = -cy(i) / dataP->MATgrid(i); 

	// ay_R, by_R, cy_R
	double rho0_R = rho0 - rho0_pi - 0.5*InnerProduct(sigq,sigq) + InnerProduct(lambda0, sigq);
	TDenseVector rho1_R = rho1 - rho1_pi + TransposeMultiply(lambda1, sigq);
	double rhov_R = rhov - rhov_pi + gamvx - 0.5;
	TDenseVector lambda0_R = lambda0 - sigq;
	TDenseMatrix lambda1_R = lambda1; 
	TDenseMatrix SIGMAlambda1_R = SIGMA*lambda1_R;
	TDenseMatrix KAPPA_rn_R = KAPPA + SIGMAlambda1_R;
	TDenseVector KAPPAtheta_rn_R = KAPPA*theta - SIGMA*lambda0_R;
	if (!YieldFacLoad_ODE3(ay_R,by_R,cy_R, dataP->TIPSgrid,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA_SIGMA,rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv*SIGMAv))
		return false; 
	for (int i=0; i<ay_R.Dimension(); i++)
		ay_R(i) = -ay_R(i)/dataP->TIPSgrid(i); 
	X = RepRowVector(dataP->TIPSgrid,Nfac,1); 
	for (int j=0; j<by_R.cols; j++)
		for (int i=0; i<by_R.rows; i++)
			by_R(i,j) = -by_R(i,j)/X(i,j); 
	by_R = Transpose(by_R);
	for (int i=0; i<cy_R.Dimension(); i++)
		cy_R(i) = -cy_R(i)/dataP->TIPSgrid(i); 

	// ay_TR, by_TR, cy_TR
	if (!YieldFacLoad_ODE3(ay_TR,by_TR,cy_TR, dataP->TIPSgrid,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA_SIGMA,rho0_R,rho1_R+rho1_L,rhov_R,KAPPAv_rn,SIGMAv*SIGMAv))
		return false; 
	for (int i=0; i<ay_TR.Dimension(); i++)
		ay_TR(i) = -ay_TR(i)/dataP->TIPSgrid(i); 
	for (int j=0; j<by_TR.cols; j++)
		for (int i=0; i<by_TR.rows; i++)
			by_TR(i,j) = -by_TR(i,j)/X(i,j); 
	by_TR = Transpose(by_TR); 
	for (int i=0; i<cy_TR.Dimension(); i++)
		cy_TR(i) = -cy_TR(i)/dataP->TIPSgrid(i); 

	// ay_L, by_L
	double KAPPA_L_rn = KAPPA_L + SIGMAlambda1_L; 
	double KAPPAtheta_L_rn=KAPPA_L*theta_L - SIGMA_L*lambda0_L; 
	if (!YieldFacLoad_ODE(ay_L,by_L, dataP->TIPSgrid, TDenseVector(1,KAPPAtheta_L_rn),TDenseMatrix(1,1,KAPPA_L_rn),TDenseMatrix(1,1,SIGMA_L*SIGMA_L), 0, TDenseVector(1,rhoL_L)))
		return false; 
	for (int i=0; i<ay_L.Dimension(); i++)
		ay_L(i) = -ay_L(i)/dataP->TIPSgrid(i); 
	for (int j=0; j<by_L.cols; j++)
		for (int i=0; i<by_L.rows; i++)
			by_L(i,j) = -by_L(i,j)/dataP->TIPSgrid(j); 
	by_L = Transpose(by_L); 

	// af, bf, cf
	TDenseVector ForecastHor(2);  
	ForecastHor(0) = 0.5; ForecastHor(1) = 1.0;
	int Nhor = ForecastHor.Dimension(); 
	af.Zeros(Nhor); 
	bf.Zeros(Nhor, Nfac); 
	cf.Zeros(Nhor); 
	for (int i=0; i<Nhor; i++)
	{
		X = MatrixExp(-ForecastHor(i)*KAPPA); 
		af(i) = ay(0) + InnerProduct(by.RowVector(1), theta, Identity(Nfac)-X) + cy(1)*(1.0-exp(-ForecastHor(i)*KAPPAv))*thetav; 
		bf.InsertRowMatrix(i,0,by.RowVector(0)*X); 
		cf(i) = cy(0)*exp(-ForecastHor(i)*KAPPAv); 
	}

	// A, B 
	TDenseMatrix Bx = MatrixExp(-dataP->Dt*KAPPA); 
	TDenseVector Ax = (Identity(Nfac)-Bx) * theta; 
	double Bv = exp(-dataP->Dt*KAPPAv); 
	double Av = (1.0-Bv) * thetav; 
	double BL = exp(-dataP->Dt*KAPPA_L); 
	double AL = (1.0-BL) * theta_L;
	B.Zeros(Nfac+3,Nfac+3); 
	B(0,0) = 1.0; 
	B.InsertRowMatrix(0,1,rho1_pi*dataP->Dt); 
	B(0,Nfac+1) = rhov_pi*dataP->Dt; 
	B.Insert(1,1,Bx); 
	B(Nfac+1,Nfac+1) = Bv; 
	B(Nfac+2,Nfac+2) = BL; 

	A.Zeros(Nfac+3); 
	A(0) = rho0_pi*dataP->Dt; 
	A.Insert(1,Ax); 
	A(Nfac+1) = Av; 
	A(Nfac+2) = AL; 

	// aI_Q, bI_Q, cI_Q 
	aI_Q.Zeros(dataP->MATgrid_options.Dimension()); 
	bI_Q.Zeros(dataP->MATgrid_options.Dimension(),Nfac); 
	cI_Q.Zeros(dataP->MATgrid_options.Dimension()); 
	double MAT, KAPPAv_Q, thetav_Q, rho0_Q, rhov_Q, tmp_aI_Q, tmp_cI_Q; 
	TDenseVector tmp_ay, tmp_cy, theta_Q, rho1_Q,  tmp_bI_Q; 
	TDenseMatrix tmp_by; 
	TDenseMatrix KAPPA_Q = KAPPA+SIGMAlambda1, KAPPA_Q_inv; 
	try {
		KAPPA_Q_inv = Inverse(KAPPA_Q);
	}
	catch(...) {
		return false; 
	}
	for (int i=0; i<dataP->MATgrid_options.Dimension(); i++)
	{
		MAT = dataP->MATgrid_options(i); 
		if (!YieldFacLoad_ODE3(tmp_ay,tmp_by,tmp_cy, TDenseVector(1,MAT),KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA_SIGMA,rho0,rho1,rhov,KAPPAv_rn,SIGMAv*SIGMAv))
			return false; 	
		theta_Q = KAPPA_Q_inv * (KAPPA*theta-SIGMA*lambda0+SIGMA_SIGMA*tmp_by.ColumnVector(0));
		KAPPAv_Q = KAPPAv + sqrt(1.0-rho*rho)*SIGMAv*gamv + rho*SIGMAv*gamvx - SIGMAv*SIGMAv*tmp_cy(0);
    		thetav_Q = KAPPAv*thetav/KAPPAv_Q;
		rho0_Q = rho0_pi - InnerProduct(sigq,lambda0) + InnerProduct(sigq, tmp_by.ColumnVector(0), Transpose(SIGMA)); 
		rho1_Q = rho1_pi - TransposeMultiply(lambda1,sigq);
    		rhov_Q = rhov_pi - gamvx + rho*SIGMAv*tmp_cy(0);
		
		if (!InfExpFacLoad_DKWv_option(tmp_aI_Q,tmp_bI_Q,tmp_cI_Q, KAPPA_Q,SIGMA,theta_Q,KAPPAv_Q,SIGMAv,thetav_Q,sigq,rho0_Q,rho1_Q,rhov_Q,rho,MAT) )
			return false; 
		aI_Q(i) = tmp_aI_Q; 
		bI_Q.InsertRowMatrix(i,0,tmp_bI_Q); 
		cI_Q(i) = tmp_cI_Q; 
	}

	// OMEGA_x
	TDenseMatrix tmp_kron; 
	try {
		tmp_kron = Inverse(Kron(-1.0*KAPPA, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA)); 
	}
	catch(...) {
		return false; 
	}
	TDenseVector tmp_x_vector(Nfac*Nfac,0.0); 
	for (int j=0; j<SIGMA_SIGMA.cols; j++)
		tmp_x_vector.Insert(j*Nfac, SIGMA_SIGMA.ColumnVector(j)); 
	tmp_x_vector = -1.0 * tmp_kron * tmp_x_vector; 
	TDenseMatrix OMEGA_x_ss(Nfac,Nfac,0.0); 
	for (int j=0; j<Nfac; j++)
		OMEGA_x_ss.InsertColumnMatrix(0,j,tmp_x_vector,j*Nfac,(j+1)*Nfac-1); 

	TDenseMatrix OMEGA_x = Bx * SIGMA_SIGMA * Transpose(Bx) - SIGMA_SIGMA; 
	tmp_x_vector.Zeros(Nfac*Nfac); 
	for (int j=0; j<OMEGA_x.cols; j++)
		tmp_x_vector.Insert(j*OMEGA_x.rows,OMEGA_x.ColumnVector(j)); 
	tmp_x_vector = tmp_kron * tmp_x_vector; 
	OMEGA_x.Zeros(Nfac,Nfac); 
	for (int j=0; j<Nfac; j++)
		OMEGA_x.InsertColumnMatrix(0,j,tmp_x_vector,j*Nfac,(j+1)*Nfac-1); 

	double OMEGA_v_ss = thetav*SIGMAv*SIGMAv/(2.0*KAPPAv);
	double OMEGA_L_ss = SIGMA_L*SIGMA_L/(2*KAPPA_L);
	double OMEGA_L = 1.0/(2.0*KAPPA_L)*(1.0-exp(-2.0*KAPPA_L*dataP->Dt))*SIGMA_L*SIGMA_L;

	double OMEGA_q_ss = 1e6;      
	TDenseMatrix KAPPA_inv; 
	try {
		KAPPA_inv = Inverse(KAPPA); 
	}	
	catch(...) {
		return false; 
	}
	TDenseVector OMEGA_xq_ss = KAPPA_inv * SIGMA*sigq; 
	double OMEGA_vq_ss = rho*SIGMAv*thetav/KAPPAv;

	TDenseMatrix OMEGA_ss(Nfac+3,Nfac+3,0.0); 
	OMEGA_ss(0,0) = OMEGA_q_ss; 
	OMEGA_ss.InsertRowMatrix(0,1,OMEGA_xq_ss); 
	OMEGA_ss(0,Nfac+1) = OMEGA_vq_ss; 
       	OMEGA_ss.InsertColumnMatrix(1,0, OMEGA_xq_ss); 
	OMEGA_ss.Insert(1,1,OMEGA_x_ss); 
	OMEGA_ss(Nfac+1,0) = OMEGA_vq_ss; 
	OMEGA_ss(Nfac+1,Nfac+1) = OMEGA_v_ss; 
	OMEGA_ss(Nfac+2,Nfac+2) = OMEGA_L_ss; 
	
	// x0_0, P0_0
	int i=0;
        while (i<dataP->p_vec.rows && dataP->p_vec(i,0) == MINUS_INFINITY)
                i++;
        x0_0.Zeros(theta.Dimension()+3);
        x0_0(0) = dataP->p_vec(i,0);
        x0_0.Insert(1,theta);
	x0_0(theta.Dimension()+1) = thetav; 
	x0_0(theta.Dimension()+2) = theta_L; 

        P0_0 = OMEGA_ss;

	// Fixed part of OMEGA
	TDenseVector OMEGA_xq = KAPPA_inv * (Identity(Nfac)-Bx) * SIGMA * sigq; 
	OMEGA.Zeros(Nfac+3,Nfac+3); 
	OMEGA.InsertRowMatrix(0,1,OMEGA_xq); 
	OMEGA.InsertColumnMatrix(1,0,OMEGA_xq); 
	OMEGA.Insert(1,1,OMEGA_x); 
	OMEGA(Nfac+2,Nfac+2) = OMEGA_L; 

	return true; 
}

bool CAR_DKWlv_o::UpdateABOmega(const TDenseVector &xtm1_tm1)
{
	double OMEGA_q = dataP->Dt*(InnerProduct(sigq,sigq) + thetav) + (xtm1_tm1(Nfac+1)-thetav)*(1-exp(-KAPPAv*dataP->Dt))/KAPPAv;
    	double OMEGA_v = thetav*SIGMAv*SIGMAv*(1.0-exp(-KAPPAv*dataP->Dt))*(1.0-exp(-KAPPAv*dataP->Dt))/(2.0*KAPPAv) + xtm1_tm1(Nfac+1)*SIGMAv*SIGMAv*(exp(-KAPPAv*dataP->Dt)-exp(-2.0*KAPPAv*dataP->Dt))/KAPPAv;
    	double OMEGA_vq = rho*SIGMAv*((1.0-exp(-KAPPAv*0.5*dataP->Dt))*(1.0-exp(-KAPPAv*0.5*dataP->Dt))/KAPPAv*thetav + xtm1_tm1(Nfac+1)*(exp(-KAPPAv*0.5*dataP->Dt)-exp(-KAPPAv*dataP->Dt))/(KAPPAv*0.5));
	
	OMEGA(0,0) = OMEGA_q; 
	OMEGA(0,Nfac+1) = OMEGA_vq;
	OMEGA(Nfac+1,0) = OMEGA_vq; 
	OMEGA(Nfac+1,Nfac+1) = OMEGA_v; 
	return true; 
}

TDenseVector CAR_DKWlv_o::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	TDenseVector model_IE_options_vector;
        if (dataP->p_vec(j,0) == MINUS_INFINITY)
	{
		b.Zeros(dataP->MATgrid.Dimension(), by.cols+3); 
		b.Insert(0,1,by); 
		b.InsertColumnMatrix(0,by.cols+1,cy); 
		a.CopyContent(ay); 
		y_obs.CopyContent(dataP->y_vec.RowVector(j));
                R.Zeros(delta_y.Dimension(), delta_y.Dimension());
                for (int i=0; i<delta_y.Dimension(); i++)
                        R(i,i) = delta_y(i) * delta_y(i);
	}
	else 
	{
		b.Zeros(dataP->MATgrid.Dimension()+1,by.cols+3); 
		b(0,0) = 1.0; 
		b.Insert(1,1,by); 
		b.InsertColumnMatrix(1,by.cols+1,cy);
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
		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+bf.rows, tmp_b.cols); 
		b.Insert(0,0,tmp_b); 
		b.Insert(tmp_b.rows, 1, bf); 
		b.InsertColumnMatrix(tmp_b.rows, bf.cols+1, cf); 

		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+af.Dimension()); 
		a.Insert(0,tmp_a); 
		a.Insert(tmp_a.Dimension(), af); 
 
		TDenseVector tmp_y_obs(y_obs);
                tmp_y_obs.CopyContent(y_obs);
                y_obs.Zeros(y_obs.Dimension()+dataP->bcf_vec.cols);
                y_obs.Insert(0,tmp_y_obs);
                y_obs.Insert(tmp_y_obs.Dimension(), dataP->bcf_vec.RowVector(j));

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
                TDenseMatrix W = 0.2* InverseMultiply(KAPPA, (MatrixExp(-hor*KAPPA)-MatrixExp(-(hor+5.0)*KAPPA)) );
		double Wv = 0.2*(exp(-KAPPAv*hor)-exp(-KAPPAv*(hor+5.0)))/KAPPAv; 
		double afLT = ay(0) + InnerProduct(by.RowVector(0),theta, Identity(Nfac)-W) + cy(0)*(1-Wv)*thetav; 
		TDenseVector bfLT = by.RowVector(0)*W; 
		double cfLT = cy(0) * Wv; 

		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+1, tmp_b.cols); 
		b.Insert(0,0,tmp_b); 
		b.InsertRowMatrix(tmp_b.rows,1,bfLT); 
		b(tmp_b.rows,1+bfLT.Dimension()) = cfLT; 

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
		b.InsertColumnMatrix(tmp_b.rows, by_TR.cols+1, cy_TR); 
		b.Insert(tmp_b.rows, by_TR.cols+2, by_L); 

		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+ay_TR.Dimension()); 
		a.Insert(0,tmp_a); 
		a.Insert(tmp_a.Dimension(), ay_TR+ay_L); 

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

	if (dataP->IE_options(j,0) != MINUS_INFINITY)
	{
		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+bI_Q.rows,tmp_b.cols); 
		b.Insert(0,0,tmp_b); 
		b.Insert(tmp_b.rows,1,bI_Q); 
		b.InsertColumnMatrix(tmp_b.rows, 1+bI_Q.cols, cI_Q); 

		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+aI_Q.Dimension()); 
		a.Insert(0,tmp_a); 
		a.Insert(tmp_a.Dimension(),aI_Q); 

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

                model_IE_options_vector = aI_Q + bI_Q * xt_tm1.SubVector(1,3) + cI_Q*xt_tm1(4); 
	}
	return model_IE_options_vector; 
}

bool CAR_DKWlv_o::SetParametersNotConsidered()
{
	if (dataP)
	{
		position_removed_from_DKW = TIndex(Nfac*Nfac + Nfac*Nfac + Nfac + 1 + Nfac + Nfac + Nfac*Nfac + dataP->y_vec.cols + dataP->bcf_vec.cols + 1 + 1 + Nfac + Nfac);
		position_removed_from_DKW_set = true; 
	}
	else 
		position_removed_from_DKW_set = false; 
	return position_removed_from_DKW_set; 
}

int CAR_DKWlv_o::NumberParameters() const
{
	int N = CAR_DKW::NumberParameters(); 
	N -= position_removed_from_DKW.size; // minus parameters of DKW that are not used 
	N += 14; // rhov, rhov_pi, KAPPAv, SIGMAv, thetav, gamv, gamvx, rho, KAPPA_L, SIGMA_L, theta_L, rhoL_L, lambda0_L, SIGMALAMBDA1_L 
	N += dataP->MATgrid_options.Dimension(); // delta_options; 
	N += Nfac; // rho1_L
	return N; 
}

TDenseVector CAR_DKWlv_o::GetParameters()
{
	if (!if_internal_parameter_set)
	{
		TDenseVector complete_parameter(NumberParameters() + position_removed_from_DKW.size,0.0); 
		complete_parameter.Insert(0, CAR_DKW::GetParameters()); 

		TIndex index_cleaned(0,NumberParameters()+position_removed_from_DKW.size-1); 
		index_cleaned.UniqSubtract(position_removed_from_DKW); 
		TDenseVector parameter_cleaned = complete_parameter.SubVector(index_cleaned); 

		int counter = CAR_DKW::NumberParameters() - position_removed_from_DKW.size; 
		parameter_cleaned(counter) = rhov; 
		counter ++; 

		parameter_cleaned(counter) = rhov_pi; 
		counter ++; 

		parameter_cleaned(counter) = KAPPAv; 
		counter ++; 

		parameter_cleaned(counter) = SIGMAv; 
		counter ++; 

		parameter_cleaned(counter) = thetav; 
		counter ++; 

		parameter_cleaned(counter) = gamv; 
		counter ++; 

		parameter_cleaned(counter) = gamvx; 
		counter ++; 

		parameter_cleaned(counter) = rho; 
		counter ++; 

		if (delta_options.Dimension() != dataP->MATgrid_options.Dimension())
			delta_options = TDenseVector(dataP->MATgrid_options.Dimension(), MINUS_INFINITY); 
		parameter_cleaned.Insert(counter, delta_options); 
		counter += delta_options.Dimension(); 

		parameter_cleaned(counter) = KAPPA_L; 
		counter ++; 

		parameter_cleaned(counter) = SIGMA_L; 
		counter ++; 
			
		parameter_cleaned(counter) = theta_L; 
		counter ++; 

		if (rho1_L.Dimension() != Nfac)
			rho1_L = TDenseVector(Nfac, MINUS_INFINITY); 
		parameter_cleaned.Insert(counter, rho1_L); 
		counter += rho1_L.Dimension(); 

		parameter_cleaned(counter) = rhoL_L; 
		counter ++; 

		parameter_cleaned(counter) = lambda0_L;
		counter ++; 
		
		parameter_cleaned(counter) = SIGMAlambda1_L; 
		counter ++; 
		
		internal_parameter.CopyContent(parameter_cleaned); 
		if_internal_parameter_set = true;
	}
	return internal_parameter; 
	
}

bool CAR_DKWlv_o::SetAsFixed(const TDenseVector &v, const std::string &which_one)
{
	TIndex local_fixed_index;
        int offset;
	if (iequals(which_one, std::string("delta_tips")) )	// because delta_tips is after sigqx, which is depleted here
	{
		offset = (CAR_DKW::NumberParameters() - position_removed_from_DKW.size) - dataP->ytips_vec.cols;
                local_fixed_index = FixedVariableParameter(delta_options, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;

	}
        else if (iequals(which_one, std::string("delta_options")) )
        {
                offset = CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8;
                local_fixed_index = FixedVariableParameter(delta_options, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else if (iequals(which_one, std::string("rho1_L")) )
        {
                offset = CAR_DKW::NumberParameters() - position_removed_from_DKW.size + dataP->MATgrid_options.Dimension() + 11; 
                local_fixed_index = FixedVariableParameter(rho1_L, v, offset);
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
        else
		return CAR_DKW::SetAsFixed(v, which_one);
}

bool CAR_DKWlv_o::SetAsFixed(double v, const std::string &which_one)
{
	TIndex local_fixed_index; 
	if (iequals(which_one, std::string("sigqx")) )
		return false; 
	else if (iequals(which_one, std::string("delta_p")) ) // because delta_p is after sigqx
	{
		delta_p = v; 
		local_fixed_index += (CAR_DKW::NumberParameters() - position_removed_from_DKW.size) - dataP->ytips_vec.cols -1; 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;	
	}
	else if (iequals(which_one, std::string("rhov")) )
	{
		rhov = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size; 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("rhov_pi")) )
	{
		rhov_pi = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 1; 
		FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("KAPPAv")) )
	{
		KAPPAv = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 2;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("SIGMAv")) )
	{
		SIGMAv = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 3;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("thetav")) )
	{
		thetav = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 4;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	} 
	else if (iequals(which_one, std::string("gamv")) )
	{
		gamv = v; 
                local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 5;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else if (iequals(which_one, std::string("gamvx")) )
	{
		gamvx = v; 
                local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 6;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else if (iequals(which_one, std::string("rho")) )
	{
		rho = v; 
                local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 7;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else if (iequals(which_one, std::string("KAPPA_L")) )
	{
		KAPPA_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension();
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else if (iequals(which_one, std::string("SIGMA_L")) )
	{
		SIGMA_L = v; 
                local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension() + 1;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
        }
	else if (iequals(which_one, std::string("theta_L")) )
	{
		theta_L = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension() + 2;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("rhoL_L")) )
	{
		rhoL_L = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension() + 2 + Nfac;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("lambda0_L")) )
	{
		lambda0_L = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension() + 2 + Nfac +1 ;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else if (iequals(which_one, std::string("SIGMAlambda1_L")) )
	{
		SIGMAlambda1_L = v; 
		local_fixed_index += CAR_DKW::NumberParameters() - position_removed_from_DKW.size + 8 + dataP->MATgrid_options.Dimension() + 2 + Nfac + 2 ;
                FixedIndex.UniqMerge(local_fixed_index);
                if_variable_index_set = false;
                return true;
	}
	else 
		return CAR_DKW::SetAsFixed(v, which_one);
}
