#include "matrix_operation.hpp"


TDenseMatrix MatrixExp(const TDenseMatrix &M)
{
	if (M.rows <=0 || M.cols <= 0 || M.rows != M.cols)
		throw dw_exception("MatrixExp() : dimension error"); 
	
	TDenseMatrix workM; 
	if (!M.column_major)
		workM = Transpose(M); 
	else 
		workM = M; 

	int ideg = 6; 
	int m = workM.cols; 
	double t = 1.0; 
	int ldh = workM.rows; 
	int lwsp = 4*m*m+ideg+1; 
	double *wsp = new double[lwsp]; 
	int *ipiv = new int[m]; 
	int iexph; 
	int ns; 
	int iflag;

	dgpadm_(&ideg, &m, &t, workM.matrix, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag); 
	
	if (iflag < 0)
		throw dw_exception("MatrixExp() : error occurred while calling gdpadm"); 

	TDenseMatrix outputM(workM.rows, workM.cols); 
	memcpy(outputM.matrix, wsp+iexph-1, sizeof(double)*workM.rows*workM.cols); 
	delete [] wsp; 
	delete [] ipiv; 
	if (!M.column_major)
		return Transpose(outputM); 
	else 
		return outputM; 
}

TDenseMatrix RepMat(const TDenseMatrix &M , int m, int n)
{
	TDenseMatrix X(M.rows*m, M.cols*n); 
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
			X.Insert(i*M.rows, j*M.cols, M); 
	}
	return X; 
}

TDenseMatrix RepRowVector(const TDenseVector &v, int m, int n)
{
	TDenseMatrix X(m, v.Dimension()*n); 
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
			X.InsertRowMatrix(i, j*v.Dimension(), v); 
	}
	return X; 
}

TDenseMatrix RepColumnVector(const TDenseVector &v, int m, int n)
{
	TDenseMatrix X(m*v.Dimension(), n); 
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
			X.InsertColumnMatrix(i*v.Dimension(), j, v); 
	}
	return X; 
}
