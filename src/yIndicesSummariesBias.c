#include <R.h>
#include "Rinternals.h"

SEXP yIndicesSummariesBias(SEXP sy, SEXP sbias, SEXP sZ, SEXP sk)
{
    int n, ldZ, i, j;
    int k = asInteger(sk);
    SEXP val;
    double *y, *bias, *ybar, *S2;
    int *N, *Z;
	
    if (k <= 0)
		error("'k' must be positive");
    if (TYPEOF(sy) != REALSXP)
		error("'y' must be of type 'double'");
    if (TYPEOF(sbias) != REALSXP)
		error("'bias' must be of type 'double'");
    if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'");
	
    n = LENGTH(sy);
    ldZ = LENGTH(sZ) / k;
	
    if (ldZ < n)
		error("'Z' has too few rows");
    if (ldZ * k != LENGTH(sZ))
		error("length of 'Z' is ot a multiple of 'k'");
	
    y = REAL(sy);
    bias = REAL(sbias);
    Z = INTEGER(sZ);
	
    PROTECT(val = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(val, 0, allocVector(INTSXP, k));
    SET_VECTOR_ELT(val, 1, allocVector(REALSXP, k));
    SET_VECTOR_ELT(val, 2, allocVector(REALSXP, k));
    N = INTEGER(VECTOR_ELT(val, 0));
    ybar = REAL(VECTOR_ELT(val, 1));
    S2 = REAL(VECTOR_ELT(val, 2));
	
    for (j = 0; j < k; j++) {
		int Nj = 0;
		double S1j = 0.0, S2j = 0.0, ybarj;
		int *z = Z + j * ldZ;
		
		for (i = 0; i < n; i++) {
			if (z[i] != 0) {
				Nj++;
				S1j = S1j + y[i] - bias[i];
			}
		}
		ybarj = Nj > 0 ? S1j / Nj : 0.0;
		N[j] = Nj;
		ybar[j] = ybarj;
		
		for (i = 0; i < n; i++) {
			if (z[i] != 0) {
				double res = y[i] - bias[i] - ybarj;
				S2j += res * res;
			}
		}
		S2[j] = S2j;
    }
	
    UNPROTECT(1);
    return val;
}
