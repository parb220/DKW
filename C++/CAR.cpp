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

double CAR::LogLikelihood(TDenseVector &logL_vec, TDenseMatrix &state, const TDenseVector &parameter)
{
	if (!SetParameters(parameter) )
		return MINUS_INFINITY;  
	return LogLikelihood(logL_vec, state); 
}

double CAR::LogLikelihood(TDenseVector &logL_vec, TDenseMatrix &state)
{
	if (dataP==NULL || dataP->NumberObservations() <=0)
		 throw dw_exception("CAR::LogLikelihood() : invalid observations");
	if (P0_0.rows != x0_0.dim || P0_0.cols != x0_0.dim)
		throw dw_exception("CAR::LogLikelihood() : dimension mismatch"); 
	
	logL_vec.Zeros(dataP->NumberObservations()); 
	state.Zeros(dataP->NumberObservations(), x0_0.dim); 

	TDenseVector xtm1_tm1 = x0_0, xt_tm1, yt_tm1; 
	TDenseMatrix Ptm1_tm1 = P0_0, Pt_tm1, Vt_tm1, Vt_tm1_inverse;   

	TDenseVector et; 
	
	double rcVt, logL=0.0; 
	for (int j=0; j<dataP->NumberObservations(); j++)
	{
		xt_tm1 = A + B*xtm1_tm1; 
		Pt_tm1 = MultiplyTranspose(B*Ptm1_tm1, B) + OMEGA; 
		
		ObservationEquation(j); 
		
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
		}
		catch (...)
		{
			logL_vec.Zeros(0); 
                        state.Zeros(0,0); 
                        return MINUS_INFINITY;
		}
		logL_vec[j] = -0.5*(LogAbsDeterminant(Vt_tm1)+InnerProduct(et, et, Vt_tm1_inverse)); 
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


void YieldFacLoad(TDenseVector &ay, TDenseMatrix &by, const TDenseMatrix &KAPPA, const TDenseMatrix &SIGMA, const TDenseVector &theta, double rho0, const TDenseVector &rho1, const TDenseVector &lambda0, const TDenseMatrix &SIGMAlambda1, const TDenseVector &Maturity)
{
	int Nfac = theta.Dimension(); 
	
	TDenseMatrix KAPPA_rn = KAPPA + SIGMAlambda1; 
	TDenseVector KAPPAtheta_rn = KAPPA*theta - SIGMA*lambda0;  

	TDenseVector A; 
	TDenseMatrix B; 
	YieldFacLoad_ODE(A, B, Maturity, KAPPAtheta_rn, KAPPA_rn, MultiplyTranspose(SIGMA, SIGMA), rho0, rho1);

	ay.Zeros(A.Dimension()); 
	for (int i=0; i<ay.Dimension(); i++)
		ay(i) = -A(i)/Maturity(i); 

	TDenseMatrix  Bx = RepRowVector(Maturity, Nfac, 1); 
	by.Zeros(B.rows, B.cols); 
	for (int j=0; j<by.cols; j++)
		for (int i=0; i<by.rows; i++)
			by(i,j) = -B(i,j)/Bx(i,j); 
	by = Transpose(by); 
}

void YieldFacLoad_ODE(TDenseVector &A, TDenseMatrix &B, const TDenseVector &TAUgrid, const TDenseVector &k, const TDenseMatrix &K, const TDenseMatrix &H, double rho0, const TDenseVector &rho1)
{
	int N = k.Dimension(); 
	int Ntau = TAUgrid.Dimension(); 

	A.Zeros(Ntau); 
	B.Zeros(N, Ntau); 

	std::vector<TDenseMatrix> expk(Ntau); 
	for (int i=0; i<Ntau; i++)
		expk[i] = MatrixExp(-TAUgrid(i) * Transpose(K)); 

	TDenseMatrix _KtInv = Inverse(-1.0*Transpose(K)); 
	TDenseVector tmp = _KtInv * (-rho1); 
	TDenseMatrix Kkron = Inverse(Kron(-1.0*K, Identity(N))+Kron(Identity(N), -1.0*K));
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
}

CAR* CAR_MinusLogLikelihood_NPSOL::model; 
void* CAR_MinusLogLikelihood_NPSOL::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
	TDenseVector logL; 
	TDenseMatrix state; 
	TDenseVector xVec(*n); 
	memcpy(xVec.vector, x, *n*sizeof(double)); 
	*f = -1.0 * model->LogLikelihood(logL, state, xVec); 
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
	double _loglikelihood_old =  LogLikelihood(logL, state, x_0), _loglikelihood_new; 
	TDenseVector delta(x_0.Dimension()), xPlusDelta; 
	for (int iter=0; iter<max_optimization_iteration; iter++)
	{
		delta.RandomNormal(); 
		for (int counter=0; counter < max_perturbation_iteration; counter ++)
		{
			xPlusDelta = x_0 + perturbation_scale * delta; 
			_loglikelihood_new = LogLikelihood(logL, state, xPlusDelta); 
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
