#include <cmath>
#include <ctgmath>
#include "prcsn.h"
#include "dw_dense_matrix.hpp"
#include "CAR.hpp"
#include "CData_FRBA.hpp"
#include "matrix_operation.hpp"

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

		try {
			rcVt = RCond(Vt_tm1); 
		}
		catch(...)
		{
			throw dw_exception("CAR::LogLikelihood(): Error occurred while evaluating the reciprocal of the condition number of Vt_tm1"); 
		}
		if (rcVt < MACHINE_EPSILON || !std::isfinite(rcVt))
		{
			logL_vec.Zeros(0); 
			state.Zeros(0,0); 
			return MINUS_INFINITY;
		} 
	
		Vt_tm1_inverse = Inverse(Vt_tm1); 
		logL_vec[j] = -0.5*(LogAbsDeterminant(Vt_tm1)+InnerProduct(et, et, Vt_tm1_inverse)); 
		logL += logL_vec[j]; 

		// update prediction of xt
		xtm1_tm1 = xt_tm1 + Pt_tm1*TransposeMultiply(b, Vt_tm1_inverse)*et; 
		Ptm1_tm1 = Pt_tm1 - Pt_tm1*TransposeMultiply(b, Vt_tm1_inverse)*b*Pt_tm1; 

		state.InsertRowMatrix(j,0,xtm1_tm1); 
	}	
	logL = -logL/(double)dataP->NumberObservations(); 
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
