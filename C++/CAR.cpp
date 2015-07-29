#include <cmath>
#include <ctgmath>
#include "prcsn.h"
#include "dw_dense_matrix.hpp"
#include "CAR.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"

extern "C" 
{
	void npoptn_(char *, int);
	void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
}

double CAR::LogLikelihood(TDenseVector &logL_vec, TDenseMatrix &state, TDenseMatrix &model_IE_options, const TDenseVector &parameter)
{
	if (!SetParameters(parameter) )
		return MINUS_INFINITY;  
	return LogLikelihood(logL_vec, state, model_IE_options); 
}

double CAR::LogLikelihood(TDenseVector &logL_vec, TDenseMatrix &state, TDenseMatrix &model_IE_options)
{
	if (dataP==NULL || dataP->NumberObservations() <=0)
		 throw dw_exception("CAR::LogLikelihood() : invalid observations");
	if (P0_0.rows != x0_0.dim || P0_0.cols != x0_0.dim)
		throw dw_exception("CAR::LogLikelihood() : dimension mismatch"); 
	
	logL_vec.Zeros(dataP->NumberObservations()); 
	state.Zeros(dataP->NumberObservations(), x0_0.dim); 
	model_IE_options.Zeros(dataP->NumberObservations(), dataP->MATgrid_options.Dimension()); 

	TDenseVector xtm1_tm1 = x0_0, xt_tm1, yt_tm1; 
	TDenseMatrix Ptm1_tm1 = P0_0, Pt_tm1, Vt_tm1, Vt_tm1_inverse;   

	TDenseVector et, model_IE_options_vector; 
	
	double logL=0.0; // rcVt; 
	for (int j=0; j<dataP->NumberObservations(); j++)
	{
		UpdateABOmega(xtm1_tm1); 
		xt_tm1 = A + B*xtm1_tm1; 
		Pt_tm1 = MultiplyTranspose(B*Ptm1_tm1, B) + OMEGA; 
		
		model_IE_options_vector = ObservationEquation(j, xt_tm1); 
		if (model_IE_options_vector.Dimension())
			model_IE_options.InsertRowMatrix(j,0,model_IE_options_vector); 
		
		yt_tm1 = a + b*xt_tm1; 
		Vt_tm1 = MultiplyTranspose(b*Pt_tm1,b) + R; 
		et = y_obs - yt_tm1; // error

		/*try {
			rcVt = RCond(Vt_tm1); 
		}
		catch(...)
		{
			throw dw_exception("CAR::LogLikelihood(): Error occurred while evaluating the reciprocal of the condition number of Vt_tm1"); 
		}
		// if (rcVt < MACHINE_EPSILON || !std::isfinite(rcVt))
		if (rcVt < 1.0e-14 || !std::isfinite(rcVt))
		// if (rcVt < 1.0e-14 || rcVt > 1.0e14)
		{
			logL_vec.Zeros(0); 
			state.Zeros(0,0); 
			return MINUS_INFINITY;
		}*/ 
	
		try {
			Vt_tm1_inverse = Inverse(Vt_tm1); 
			TDenseMatrix T = Cholesky(Vt_tm1); 
			logL_vec[j] = -0.5*(LogAbsDeterminant_Cholesky(T)+InnerProduct(et, et, Vt_tm1_inverse)); 
		}
		catch (...)
		{
			logL_vec.Zeros(0); 
                        state.Zeros(0,0); 
			model_IE_options.Zeros(0,0); 
                        return MINUS_INFINITY;
		}
		logL += logL_vec[j]; 

		// update prediction of xt
		xtm1_tm1 = xt_tm1 + Pt_tm1*TransposeMultiply(b, Vt_tm1_inverse)*et; 
		Ptm1_tm1 = Pt_tm1 - Pt_tm1*TransposeMultiply(b, Vt_tm1_inverse)*b*Pt_tm1; 

		state.InsertRowMatrix(j,0,xtm1_tm1); 
	}	
	logL = logL/(double)dataP->NumberObservations(); 
	if (!std::isfinite(logL)|| std::isnan(logL))
		return MINUS_INFINITY; 
	return logL; 
}

CAR::CAR(int _Nfac, CData_FRBA* _dataP) : 
Nfac(_Nfac), dataP(_dataP), 
A(), B(), OMEGA(),
x0_0(), P0_0(),
a(), b(), y_obs(), R()
{}

CAR::CAR(const CAR &right) :
Nfac(right.Nfac), dataP(right.dataP) 
{
	A.CopyContent(right.A); 
	B.CopyContent(right.B); 
	OMEGA.CopyContent(right.OMEGA); 
	x0_0.CopyContent(right.x0_0); 
	P0_0.CopyContent(right.P0_0); 
	a.CopyContent(right.a); 
	b.CopyContent(right.b); 
	y_obs.CopyContent(right.y_obs); 
	R.CopyContent(right.R); 
}


bool YieldFacLoad(TDenseVector &ay, TDenseMatrix &by, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, double rho0, const TDenseVector &rho1, const TDenseVector &lambda0, const TDenseMatrix &SIGMAlambda1, const TDenseVector &Maturity)
{
	int Nfac = theta.Dimension(); 
	
	TDenseMatrix KAPPA_rn = KAPPA + SIGMAlambda1; 
	TDenseVector KAPPAtheta_rn = KAPPA*theta - SIGMA*lambda0;  

	TDenseVector A; 
	TDenseMatrix B; 
	if(!YieldFacLoad_ODE(A, B, Maturity, KAPPAtheta_rn, KAPPA_rn, MultiplyTranspose(SIGMA, SIGMA), rho0, rho1))
		return false; 

	ay.Zeros(A.Dimension()); 
	for (int i=0; i<ay.Dimension(); i++)
		ay(i) = -A(i)/Maturity(i); 

	TDenseMatrix  Bx = RepRowVector(Maturity, Nfac, 1); 
	by.Zeros(B.rows, B.cols); 
	for (int j=0; j<by.cols; j++)
		for (int i=0; i<by.rows; i++)
			by(i,j) = -B(i,j)/Bx(i,j); 
	by = Transpose(by); 
	return true; 
}

bool YieldFacLoad_ODE(TDenseVector &A, TDenseMatrix &B, const TDenseVector &TAUgrid, const TDenseVector &k, const TDenseMatrix &K, const TDenseMatrix &H, double rho0, const TDenseVector &rho1)
{
	int N = k.Dimension(); 
	int Ntau = TAUgrid.Dimension(); 

	A.Zeros(Ntau); 
	B.Zeros(N, Ntau); 

	std::vector<TDenseMatrix> expk(Ntau); 
	for (int i=0; i<Ntau; i++)
		expk[i] = MatrixExp(-TAUgrid(i) * Transpose(K)); 

	TDenseMatrix _KtInv; 
	try {
		_KtInv = Inverse(-1.0*Transpose(K)); 
	}
	catch(...) {
		return false; 
	}
	TDenseVector tmp = _KtInv * (-rho1); 
	TDenseMatrix Kkron; 
	try { 
		Kkron = Inverse(Kron(-1.0*K, Identity(N))+Kron(Identity(N), -1.0*K));
	}
	catch(...) {
		return false; 
	}
	for (int i=0; i<Ntau; i++)
	{
		B.InsertColumnMatrix(0,i, (expk[i]-Identity(N))*tmp); 
		double tau = TAUgrid(i); 
		TDenseVector K0 = ((expk[i]-Identity(N))*_KtInv-tau*Identity(N))*tmp; 

		TDenseMatrix K2 = TransposeMultiply(expk[i], H*expk[i]) - H;  
		TDenseVector K2_vector(N*N); 
		for (int j=0; j<N; j++)
			K2_vector.Insert(j*N, K2.ColumnVector(j)); 
		K2_vector = Kkron * K2_vector; 
		K2.Zeros(N, N); 
		for (int j=0; j<N; j++)
			K2.InsertColumnMatrix(0,j, K2_vector, j*N, (j+1)*N-1); 

		TDenseMatrix K1 = H * _KtInv * (expk[i]-Identity(N)); 
		A(i) = InnerProduct(k, K0) + 0.5*InnerProduct(tmp, tmp, K2-K1-Transpose(K1)+tau*H) - tau*rho0; 
	} 
	return true; 
}

bool YieldFacLoad_ODE3(TDenseVector &A, TDenseMatrix &B, TDenseVector &C, const TDenseVector &TAUgrid, const TDenseVector &k, const TDenseMatrix &K, double h, const TDenseMatrix &H, double rho0, const TDenseVector &rho1, double rhov, double Kv, double Hv)
{
	int Ntau = TAUgrid.Dimension(); 

	TDenseVector A1; 
	if (!YieldFacLoad_ODE(A1,B,TAUgrid,k,K,H,rho0,rho1))
		return false; 

	double a = -rhov, b = -Kv, c = 0.5*Hv, ETA; 
	if (b*b - 4.0*a*c > 0)
		ETA = sqrt(b*b - 4.0*a*c); 
	else 
		return false; 
	
	C.Zeros(Ntau); 
	for (int i=0; i<C.Dimension(); i++)
		C(i) = 2.0*a*(1.0-exp(-ETA*TAUgrid(i)))/((ETA+b)*exp(-ETA*TAUgrid(i))+(ETA-b));

	TDenseVector A2(Ntau,0.0); 
	if (a != 0)
	{
		for (int i=0; i<A2.Dimension(); i++)
			A2(i) = h*2.0*a*(TAUgrid(i)/(ETA-b)+2.0/(ETA*ETA-b*b)*log(((ETA+b)*exp(-ETA*TAUgrid(i))+(ETA-b))/(2.0*ETA)));
	}
	A = A1+A2; 
	return true; 	
}

CAR* CAR_MinusLogLikelihood_NPSOL::model; 
void* CAR_MinusLogLikelihood_NPSOL::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
	TDenseVector logL; 
	TDenseMatrix state;
	TDenseMatrix model_IE_options;  
	TDenseVector xVec(*n); 
	memcpy(xVec.vector, x, *n*sizeof(double)); 
	*f = -1.0 * model->LogLikelihood(logL, state, model_IE_options, xVec); 
}

double CAR::MaximizeLogLikelihood(TDenseVector &x_0, const TDenseVector &lower_x, const TDenseVector &upper_x, const std::vector<std::string> &option, int max_perturbation_iteration, double perturbation_scale, int max_optimization_iteration)
{
	// set npoptn_
	const std::string COLD_START = std::string("Cold Start");
        const std::string NO_PRINT_OUT = std::string("Major print level = 0");
        const std::string DERIVATIVE_LEVEL = std::string("Derivative level = 0"); 
	
	CAR_MinusLogLikelihood_NPSOL::model = this; 
	
	int n = x_0.Dimension(); // number of variables 
	double *bl = new double[n], *bu = new double[n]; 
	if (lower_x.Dimension() == n)
		memcpy(bl, lower_x.vector, sizeof(double)*n); 
	else 
	{
		for (int i=0; i<n; i++)
			bl[i] = MINUS_INFINITY; 
	}

	if (upper_x.Dimension() == n)
		memcpy(bu, upper_x.vector, sizeof(double)*n); 
	else 
	{
		for (int i=0; i<n; i++)
			bu[i] = PLUS_INFINITY; 
	}

	// Linear constraint
	int nclin = 0; // number of linear constraints
	int ldA = 1; // number of rows of A matrix
	double *A = new double[ldA*1];  

	// Nonlinear constraint
	int ncnln = 0; // number of nonlinear constraints
	int ldJ = 1; // number of rows oof cJac
	double *c = new double[1], *cJac = new double[ldJ*1]; 	

	// R provides the upper-triangular Cholesky factorization of the initial approximation
	// of the Hessian of the Lagrangian
	int ldR = n; // number of rows of R
	double *R = new double[ldR*n]; 

	int *istate = new int[n+nclin+ncnln]; 
	double *clamda = new double[n+nclin+ncnln]; 
	int leniw = 3*n + nclin + 2*ncnln; 
	int *iw = new int[leniw];
	int lenw = 20 *n;  // because nclin==0 && ncnln ==0
	double *w = new double[lenw]; 

	double *x = new double[n]; 

	int inform, iter; 
	double f;
	double *g = new double[n]; 
	
	// find an appropriate initial value
	TDenseVector logL; 
	TDenseMatrix state;
	TDenseMatrix model_IE_options;  
	double _loglikelihood_old =  LogLikelihood(logL, state, model_IE_options, x_0), _loglikelihood_new; 
	TDenseVector delta(x_0.Dimension()), xPlusDelta; 
	for (int iteration=0; iteration<max_optimization_iteration; iteration++)
	{
		delta.RandomNormal(); 
		for (int counter=0; counter < max_perturbation_iteration; counter ++)
		{
			xPlusDelta = x_0 + perturbation_scale * delta; 
			_loglikelihood_new = LogLikelihood(logL, state, model_IE_options, xPlusDelta); 
			if (_loglikelihood_new > _loglikelihood_old)
			{
				x_0 = xPlusDelta; 
				_loglikelihood_old = _loglikelihood_new; 
			}
			delta = 0.5*delta; 
		}

		npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
       	 	npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
        	npoptn_((char*)COLD_START.c_str(), COLD_START.length());

		for (int i=0; i<(int)option.size(); i++)
			npoptn_((char*)option[i].c_str(), option[i].length());

		memcpy(x, x_0.vector, sizeof(double)*n);
	
		npsol_(&n,&nclin,&ncnln,&ldA,&ldJ,&ldR,A,bl,bu,NULL,CAR_MinusLogLikelihood_NPSOL::function, &inform, &iter, istate, c, cJac, clamda, &f, g, R, x, iw, &leniw, w, &lenw); 
	
		_loglikelihood_new = -f; 
		if (_loglikelihood_new > _loglikelihood_old)
		{
			memcpy(x_0.vector, x, sizeof(double)*n); 
			_loglikelihood_old = _loglikelihood_new; 
		}
	}
	
	delete [] bl; 
	delete [] bu; 
	delete [] A; 
	delete [] c; 
	delete [] cJac; 
	delete [] R; 
	delete [] istate; 
	delete [] clamda; 
	delete [] iw; 
	delete [] w; 
	delete [] x; 
	delete [] g; 

	return _loglikelihood_old; 
}

bool InfExpFacLoad(double &A, TDenseVector &B, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, const TDenseVector &sigq, double sigqx, double rho0_pi, const TDenseVector &rho1_pi, double Maturity)
{
	int Nfac = theta.Dimension(); 
	TDenseMatrix temp_kron; 
	try {
		temp_kron = Inverse(Kron(-1.0*KAPPA, Identity(Nfac))+Kron(Identity(Nfac), -1.0*KAPPA)); 
	}
	catch(...) {
		return false; 
	}

	TDenseMatrix temp_km = MatrixExp(-Maturity*KAPPA); 
	TDenseMatrix temp_x = temp_km * SIGMA; 
	TDenseMatrix before_reshape = MultiplyTranspose(temp_x, temp_x) - MultiplyTranspose(SIGMA, SIGMA); 
	TDenseVector after_reshape(Nfac*Nfac,0.0); 
	for (int j=0; j<before_reshape.cols; j++)
		after_reshape.Insert(j*before_reshape.rows, before_reshape.ColumnVector(j)); 
	after_reshape = temp_kron * after_reshape; 
	TDenseMatrix OMEGA_x(Nfac, Nfac, 0.0); 
	for (int j=0; j<Nfac; j++)
		OMEGA_x.InsertColumnMatrix(0,j,after_reshape, j*Nfac, (j+1)*Nfac-1); 

	TDenseMatrix _km_inv; 
	try {
		_km_inv = Inverse(KAPPA*Maturity); 
	}
	catch(...) {
		return false; 
	}
 	TDenseMatrix bx = _km_inv * (Identity(Nfac)-temp_km); 
	TDenseVector ax = (Identity(Nfac)-bx) * theta; 

	TDenseMatrix _k_inv; 
	try { 
		_k_inv = Inverse(KAPPA); 
	}
	catch (...)
	{
		return false; 
	}

	TDenseVector sigq2 = sigq + TransposeMultiply(_k_inv*SIGMA, rho1_pi); 
	double tempIU_1 = InnerProduct(sigq2,sigq2) + sigqx*sigqx; 
	double tempIU_2 =  - 2.0*InnerProduct(sigq2, TransposeMultiply(_km_inv*SIGMA, TransposeMultiply(_k_inv*(Identity(Nfac)-temp_km), rho1_pi) ) ); 
	double tempIU_3 = InnerProduct(rho1_pi, _k_inv*(1.0/Maturity)*OMEGA_x*TransposeMultiply(_k_inv, rho1_pi));

	double tempIU = tempIU_1 + tempIU_2 + tempIU_3; 
	B = TransposeMultiply(bx, rho1_pi); 
	A = (rho0_pi + InnerProduct(rho1_pi, ax)) + 0.5*tempIU;  
	return true; 
}

bool InfExpFacLoad_DKWv_option(double &aI, TDenseVector &bI, double &cI, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, double KAPPAv, double SIGMAv, double thetav, const TDenseVector &sigq, double rho0_pi, const TDenseVector &rho1_pi, double rhov_pi, double rho, double Maturity)
{
	int Nfac = theta.Dimension();
	TDenseMatrix tmp_kron; 
	try {
		tmp_kron = Inverse(Kron(-1.0*KAPPA,Identity(Nfac))+Kron(Identity(Nfac),-1.0*KAPPA)); 
	} 	
	catch(...) {
		return false; 
	}
	TDenseMatrix SIGMA_SIGMA = MultiplyTranspose(SIGMA,SIGMA); 
	TDenseMatrix tmp_expm = MatrixExp(-Maturity*KAPPA); 
	TDenseMatrix tmp_X = tmp_expm*SIGMA_SIGMA*Transpose(tmp_expm) - SIGMA_SIGMA; 
	TDenseVector tmp_X_vector(Nfac*Nfac,0.0); 
	for (int j=0; j<tmp_X.cols; j++)
		tmp_X_vector.Insert(j*tmp_X.rows, tmp_X.ColumnVector(j)); 
	tmp_X_vector = tmp_kron * tmp_X_vector; 
	TDenseMatrix OMEGA_x(Nfac,Nfac); 
	for (int j=0; j<Nfac; j++)
		OMEGA_x.InsertColumnMatrix(0,j,tmp_X_vector,j*Nfac,(j+1)*Nfac-1); 

	TDenseMatrix KAPPA_maturity_inv; 
	try {
		KAPPA_maturity_inv = Inverse(Maturity*KAPPA); 
	}
	catch(...) {
		return false; 
	}
	TDenseMatrix bx = KAPPA_maturity_inv * (Identity(Nfac)-tmp_expm); 
	TDenseVector ax = (Identity(Nfac)-bx) * theta; 
	double bv = 1.0/(KAPPAv*Maturity) * (1.0-exp(-KAPPAv*Maturity)); 
	double av = (1.0-bv) * thetav; 

	TDenseMatrix KAPPA_inv; 
	try {
		KAPPA_inv = Inverse(KAPPA); 
	}
	catch(...) {
		return false; 
	}
	TDenseVector sigq2 = sigq + TransposeMultiply(KAPPA_inv*SIGMA, rho1_pi); 
	double H0 = InnerProduct(sigq2,sigq2) - 2.0*InnerProduct(sigq2, rho1_pi, Transpose(KAPPA_inv*(Identity(Nfac)-tmp_expm)*KAPPA_maturity_inv*SIGMA) ) + InnerProduct(rho1_pi,rho1_pi, KAPPA_inv*(1.0/Maturity)*OMEGA_x*Transpose(KAPPA_inv)); 

	double tmp1 = rhov_pi*SIGMAv/KAPPAv; 
	double H1 = (1.0+2.0*rho*tmp1+tmp1*tmp1)-2.0*(rho+tmp1)*tmp1/(KAPPAv*Maturity)*(1.0-exp(-KAPPAv*Maturity))+tmp1*tmp1/(2.0*KAPPAv*Maturity)*(1.0-exp(-2*KAPPAv*Maturity));
	double H2 = exp(-KAPPAv*Maturity)*((1.0+2.0*rho*tmp1+tmp1*tmp1)/(KAPPAv*Maturity)*(exp(KAPPAv*Maturity)-1.0)-2.0*(rho+tmp1)*tmp1+tmp1*tmp1/(KAPPAv*Maturity)*(1.0-exp(-KAPPAv*Maturity)));
	
	bI = TransposeMultiply(bx,rho1_pi); 
	cI = bv*rhov_pi+0.5*H2; 
	aI = (rho0_pi+InnerProduct(rho1_pi,ax)+rhov_pi*av)+0.5*(H0+thetav*(H1-H2));

	return true; 
}

bool iequals(const std::string& a, const std::string& b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

TIndex FixedVariableParameter(TDenseMatrix &destination, const TDenseMatrix & source, int offset)
{
        TIndex fixed_index; 
        destination = source; 
        for (int j=0; j<destination.cols; j++)
        {
                for (int i=0; i<destination.rows; i++)
                {
                        if (destination(i,j) > MINUS_INFINITY && destination(i,j) < PLUS_INFINITY)
                                fixed_index += offset + j*destination.rows + i;
                }
        }
        return fixed_index;
}

TIndex FixedVariableParameter(TDenseVector &destination, const TDenseVector &source, int offset)
{
        TIndex fixed_index;
        destination = source;   
        for (int j=0; j<destination.Dimension(); j++)
        {
                if (destination(j) > MINUS_INFINITY && destination(j) < PLUS_INFINITY)
                        fixed_index += offset + j;
        }
        return fixed_index;
}



