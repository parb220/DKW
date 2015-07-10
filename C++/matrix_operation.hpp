#ifndef _MATRIX_OPERATON_HEADER
#define _MATRIX_OPERATION_HEADER
extern "C"
{
        void dgpadm_(int *ideg, int *m, double *t, double *H, int *ldh, double *wsp, int *lwsp, int *ipiv, int *iexph, int *ns, int *iflag);
}

#include "dw_dense_matrix.hpp"
TDenseMatrix MatrixExp(const TDenseMatrix &M); 
TDenseMatrix RepMat(const TDenseMatrix &M , int m, int n); 
TDenseMatrix RepRowVector(const TDenseVector &v, int m, int n); 
TDenseMatrix RepColumnVector(const TDenseVector &v, int m, int n); 

#endif
