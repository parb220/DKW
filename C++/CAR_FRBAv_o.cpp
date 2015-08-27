#include "CAR_FRBAv_o.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"
#include "prcsn.h"
#include "dw_rand.h"

CAR_FRBAv_o::CAR_FRBAv_o(int _Nfac, CData_FRBA *_dataP) : CAR_DKWv_o(_Nfac, _dataP), 
delta_option_position_set(false),
delta_option_index_start(), delta_option_index_end()
{
	if (dataP)
	{
		delta_option_index_start = CAR_DKW_o::NumberParameters(); 
		delta_option_index_end = CAR_DKW_o::NumberParameters()+dataP->MATgrid_options.Dimension()-1; 
		delta_option_position_set = true; 
	}
}

CAR_FRBAv_o::CAR_FRBAv_o(const CAR_FRBAv_o &right) : CAR_DKWv_o(right),
delta_option_position_set(right.delta_option_position_set),
delta_option_index_start(right.delta_option_index_start), 
delta_option_index_end(right.delta_option_index_end)
{}

bool CAR_FRBAv_o::SetDeltaOptionPosition()
{
	if (!delta_option_position_set)
	{	
		delta_option_index_start = CAR_DKW_o::NumberParameters(); 
		delta_option_index_end = CAR_DKW_o::NumberParameters()+dataP->MATgrid_options.Dimension()-1; 
		delta_option_position_set = true; 
	}
	return delta_option_position_set; 
}

int CAR_FRBAv_o::NumberParameters() const
{
	int _N = CAR_DKWv_o::NumberParameters(); 
	_N = _N - dataP->MATgrid_options.Dimension() + 1; // delta_options dealt as a scalar rather than a vector
	return _N; 
}

TDenseVector CAR_FRBAv_o::GetParameters() 
{
	if (!if_internal_parameter_set)
	{
		TDenseVector DKWv_o_parameter = CAR_DKWv_o::GetParameters(); 
		TDenseVector FRBAv_o_parameter(this->NumberParameters(),0.0); 
		
		SetDeltaOptionPosition(); 
		FRBAv_o_parameter.Insert(0, DKWv_o_parameter, 0, delta_option_index_start); //up to delta_options(0);
		FRBAv_o_parameter.Insert(delta_option_index_start+1, DKWv_o_parameter, delta_option_index_end+1, CAR_DKWv_o::NumberParameters()-1); // skipped delta_options(1) ... delta_options(end)

        	internal_parameter.CopyContent(FRBAv_o_parameter);
        	if_internal_parameter_set = true;
	}
	return internal_parameter;  
}

bool CAR_FRBAv_o::SetParameters_FromVectorToMatrix(const TDenseVector &_parameter)
{
	if (_parameter.Dimension() == CAR_DKWv_o::NumberParameters())
		return CAR_DKWv_o::SetParameters_FromVectorToMatrix(_parameter); 

	else if (_parameter.Dimension() == this->NumberParameters())
	{
		SetDeltaOptionPosition();
		TDenseVector augmented_parameter(CAR_DKWv_o::NumberParameters(),0.0);
		augmented_parameter.Insert(0, _parameter, 0, delta_option_index_start-1); 
		augmented_parameter.Insert(delta_option_index_start, TDenseVector(dataP->MATgrid_options.Dimension(), _parameter(delta_option_index_start)) ); 
		augmented_parameter.Insert(delta_option_index_end+1, _parameter, delta_option_index_start+1, _parameter.Dimension()-1); 
		return CAR_DKWv_o::SetParameters_FromVectorToMatrix(augmented_parameter); 
	}
	else 
		return false; 
}

bool CAR_FRBAv_o::SetAsFixed(double v, const std::string &which_one)
{
	if (iequals(which_one, std::string("delta_options")) )
		return SetAsFixed(TDenseVector(dataP->MATgrid_options.Dimension(), v), which_one); 		
	else 
		return CAR_DKWv_o::SetAsFixed(v, which_one); 
}

TDenseVector CAR_FRBAv_o::ObservationEquation(int j, const TDenseVector &xt_tm1)
{
	if (dataP->p_vec(j,0) == MINUS_INFINITY)
        {
                b.Zeros(dataP->MATgrid.Dimension(), by.cols+2);
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
                b.Zeros(dataP->MATgrid.Dimension()+1,by.cols+2);
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
                TDenseMatrix W = 0.2* Multiply(KAPPA_inverse, (MatrixExp(-hor*KAPPA)-MatrixExp(-(hor+5.0)*KAPPA)) );
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
                b.Zeros(tmp_b.rows+by_R.rows, tmp_b.cols);
                b.Insert(0,0,tmp_b);
                b.Insert(tmp_b.rows,1,by_R);
                b.InsertColumnMatrix(tmp_b.rows, by_R.cols+1, cy_R);

                TDenseVector tmp_a(a);
                tmp_a.CopyContent(a);
                a.Zeros(tmp_a.Dimension()+ay_R.Dimension());
                a.Insert(0,tmp_a);
                a.Insert(tmp_a.Dimension(), ay_R);

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

	if (dataP->caps_vec(j,0) != MINUS_INFINITY)
	{
		TDenseVector ay_caps, cy_caps, model_prc_caps; 
		TDenseMatrix by_caps; 
		if (!CapsFacLoad_FRBAv_o(ay_caps, by_caps, cy_caps, model_prc_caps, xt_tm1.SubVector(1,Nfac), xt_tm1(Nfac+1), dataP->MATgrid_caps, dataP->STRIKEgrid_caps, false) )
			throw dw_exception("Error in loading caps_vec"); 
		
		TDenseVector scale(6); 
		scale(0)=1.0; scale(1)=1.0; scale(2)=1.0; 
		scale(3)=3.0; scale(4)=3.0; scale(5)=3.0; 	

		TDenseVector caps_prc(dataP->caps_vec.RowVector(j)); 
		caps_prc.CopyContent(dataP->caps_vec.RowVector(j)); 
		for (int i=0; i<scale.Dimension(); i++)
		{
			ay_caps(i) = ay_caps(i)/scale(i); 
			cy_caps(i) = cy_caps(i)/scale(i); 
			caps_prc(i) = caps_prc(i)/scale(i); 	
			for (int j=0; j<by_caps.cols; j++)
				by_caps(i,j) = by_caps(i,j)/scale(i); 
		}

		TDenseVector tmp_a(a); 
		tmp_a.CopyContent(a); 
		a.Zeros(tmp_a.Dimension()+ay_caps.Dimension()); 
		a.Insert(0,tmp_a); 
		a.Insert(tmp_a.Dimension(), ay_caps); 

		TDenseMatrix tmp_b(b); 
		tmp_b.CopyContent(b); 
		b.Zeros(tmp_b.rows+by_caps.rows, tmp_b.cols); 
		b.Insert(0,0,tmp_b); 
		b.Insert(tmp_b.rows,1,by_caps); 
		b.InsertRowMatrix(tmp_b.rows,1+by_caps.cols, cy_caps); 

		TDenseVector tmp_y_obs(y_obs); 
		tmp_y_obs.CopyContent(y_obs); 
		y_obs.Zeros(tmp_y_obs.Dimension()+caps_prc.Dimension()); 
		y_obs.Insert(0,tmp_y_obs); 
		y_obs.Insert(tmp_y_obs.Dimension(), caps_prc); 

		TDenseMatrix tmp_R(R); 
		tmp_R.CopyContent(R); 
		R.Zeros(tmp_R.rows+dataP->MATgrid_caps.Dimension(), tmp_R.cols+dataP->MATgrid_caps.Dimension()); 
		R.Insert(0,0,tmp_R); 
		for (int i=0; i<dataP->MATgrid_caps.Dimension(); i++)
			R(i+tmp_R.rows, i+tmp_R.cols) = delta_options(0)*delta_options(0); 	
	}

	return TDenseVector(); 	
}

bool CAR_FRBAv_o::CapsFacLoad_FRBAv_o(TDenseVector &ay_caps, TDenseMatrix &by_caps, TDenseVector &cy_caps, TDenseVector &model_prc_caps, const TDenseVector &xt_tm1, double vt_tm1, const TDenseVector &Maturity, const TDenseVector &Strike, bool floor_flag)
{
	int Nmat = Maturity.Dimension(); 
	int Nk = Strike.Dimension(); 
	int Ncaps = Nmat * Nk; 

	// ay, by, cy already calculated 
	ay_caps.Zeros(Ncaps); 
	by_caps.Zeros(Ncaps,Nfac); 
	cy_caps.Zeros(Ncaps); 
	model_prc_caps.Zeros(Ncaps); 

	TDenseMatrix lambda1 = SIGMA_inverse * SIGMAlambda1; 
    	TDenseVector rho1_Q = rho1_pi - TransposeMultiply(lambda1,sigq);
	TDenseVector sigq2 = sigq + TransposeMultiply(Inv_KAPPA_rn*SIGMA, rho1_Q);
	for (int i_tau=0; i_tau<Nmat; i_tau++)
	{
		double tau = Maturity(i_tau); 
		// KAPPA_Q = KAPPA_rn
		TDenseVector  theta_Q = Inv_KAPPA_rn * ((KAPPA*theta - SIGMA*lambda0 - MultiplyTranspose(SIGMA,SIGMA)*by.RowVector(i_tau))*tau);
		
		double KAPPAv_Q = KAPPAv + sqrt(1.0-rho*rho)*SIGMAv*gamv + rho*SIGMAv*gamvx + SIGMAv*SIGMAv*cy(i_tau)*tau;
     		double thetav_Q = KAPPAv*thetav/KAPPAv_Q; 

		double rho0_Q = rho0_pi - InnerProduct(sigq,lambda0) - InnerProduct(sigq, by.RowVector(i_tau), Transpose(SIGMA))*tau;
    		double rhov_Q = rhov_pi - gamvx - rho*SIGMAv*cy(i_tau)*tau;

		TDenseMatrix Bx = MatrixExp(-1.0*KAPPA_rn*tau); 
		TDenseMatrix bx_Q = (1.0/tau)*Inv_KAPPA_rn*(Identity(Nfac)-Bx);
		TDenseVector ax_Q = (Identity(Nfac)-bx_Q)*theta_Q;	
    		double bv_Q = 1.0/(KAPPAv_Q*tau)*(1.0-exp(-KAPPAv_Q*tau));
		double av_Q = (1.0-bv_Q)*thetav_Q;

		double aI_QQ = rho0_Q + InnerProduct(rho1_Q,ax_Q)+ rhov_Q*av_Q;
    		TDenseVector bI_QQ = TransposeMultiply(bx_Q, rho1_Q); 
    		double cI_QQ = rhov_Q*bv_Q;	

		TDenseMatrix Xx = Bx*SIGMA; 
		Xx = MultiplyTranspose(Xx, Xx) - MultiplyTranspose(SIGMA, SIGMA); 
		TDenseVector Xx_vector(Nfac*Nfac,0.0); 
		for (int j=0; j<Xx.cols; j++)
			Xx_vector.Insert(j*Xx.rows, Xx.ColumnVector(j)); 
		Xx_vector = Inv_Kron_KAPPA_rn * Xx_vector; 	
		TDenseMatrix OMEGA_x_Q(Nfac,Nfac,0.0); 
		for (int i=0; i<Nfac; i++)
			OMEGA_x_Q.InsertColumnMatrix(0,i,Xx_vector, i*Nfac,(i+1)*Nfac-1); 

		double H0 = InnerProduct(sigq2,sigq2) - 2.0*InnerProduct(Inv_KAPPA_rn* ((Identity(Nfac)-Bx)*(1.0/tau)*(Inv_KAPPA_rn*(SIGMA*sigq2))), rho1_Q)  + InnerProduct(rho1_Q, Inv_KAPPA_rn*(((1.0/tau)*OMEGA_x_Q)* (Transpose(Inv_KAPPA_rn)*rho1_Q)) ); 
		double tmp1 = rhov_Q*SIGMAv/KAPPAv_Q;
    		double H1 = (1.0 + 2.0*rho*tmp1+tmp1*tmp1) - 2.0*(rho+tmp1)*tmp1/(KAPPAv_Q*tau)*(1-exp(-1.0*KAPPAv_Q*tau)) + tmp1*tmp1/(2.0*KAPPAv_Q*tau)*(1.0-exp(-2.0*KAPPAv_Q*tau));
    		double H2 = exp(-1.0*KAPPAv_Q*tau)*((1.0+2.0*rho*tmp1+tmp1*tmp1)/(KAPPAv_Q*tau)*(exp(KAPPAv_Q*tau)-1.0) - 2.0*(rho+tmp1)*tmp1+tmp1*tmp1/(KAPPAv_Q*tau)*(1.0-exp(-KAPPAv_Q*tau)));

    		double dI_Q = H0 + thetav_Q * (H1-H2);
    		double eI_Q = H2;
	
		for (int i_k=0; i_k<Nk; i_k++)
		{
			double strike_CAPS = Strike(i_k); 
			if (strike_CAPS <= -1.0)
				return false; 

			double tmp_mean = aI_QQ + InnerProduct(bI_QQ, xt_tm1) + cI_QQ * vt_tm1; 
			double tmp_var = dI_Q + eI_Q*vt_tm1;
			double h0 = tau*(tmp_mean+tmp_var*0.5);
        		double h1 = (-log(1.0+strike_CAPS)+(tmp_mean+tmp_var))/sqrt(tmp_var/tau);
        		double h2 = (-log(1.0+strike_CAPS)+tmp_mean)/sqrt(tmp_var/tau);

        		TDenseVector h0x = tau*bI_QQ; 
			TDenseVector h1x = bI_QQ*(1.0/sqrt(tmp_var/tau)); 
			TDenseVector h2x = h1x;
        		double h0v = tau*(cI_QQ+eI_Q/2.0);
        		double h1v = (cI_QQ+eI_Q)/sqrt(tmp_var/tau) - 0.5*eI_Q/tmp_var*h1;
        		double h2v = cI_QQ/sqrt(tmp_var/tau)-0.5*eI_Q/tmp_var*h2;

			if (floor_flag) 
			{
				h1 = -h1; 
				h1x = -1.0*h1x; 
				h1v = -h1v; 

				h2 = -h2; 
				h2x = -1.0*h2x; 
				h2v = -h2v; 
			}

			double tmp_a = log(exp(h0)*dw_normal_cdf(h1)- pow(1.0+strike_CAPS, tau)*dw_normal_cdf(h2));
        		TDenseVector tmp_b = (exp(h0)*(dw_normal_cdf(h1)*h0x + dw_normal_pdf(h1)*h1x)-pow(1.0+strike_CAPS, tau)*dw_normal_pdf(h2)*h2x)*(1.0/exp(tmp_a));
        		double tmp_c = (exp(h0)*(dw_normal_cdf(h1)*h0v + dw_normal_pdf(h1)*h1v)-pow(1+strike_CAPS,tau)*dw_normal_pdf(h2)*h2v)/exp(tmp_a);

        		int id = i_k+(i_tau-1)*Nk;
        		ay_caps(id) = -tau*ay(i_tau) + tmp_a- InnerProduct(tmp_b, xt_tm1)- tmp_c*vt_tm1;
        		by_caps.InsertRowMatrix(id,0, -tau*by.RowVector(i_tau)+tmp_b);
        		cy_caps(id) = -tau*cy(i_tau) + tmp_c;

        		model_prc_caps(id) = -tau*(ay(i_tau)+InnerProduct(by.RowVector(i_tau),xt_tm1)+cy(i_tau)*vt_tm1)+tmp_a;
		}	
	}	
	return true; 
}
