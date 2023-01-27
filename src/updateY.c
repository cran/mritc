#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP updateY(SEXP ss, SEXP ssubvox, SEXP sZ, SEXP smu, SEXP ssigma)
{
	if (TYPEOF(ss) != REALSXP)
		error("'s' must be of type 'double'.");
    if (TYPEOF(ssubvox) != INTSXP)
		error("'subvox' must be of type 'integer'.");
	if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'.");
	if (TYPEOF(smu) != REALSXP)
		error("'mu' must be of type 'double'.");
    if (TYPEOF(ssigma) != REALSXP)
		error("'sigma' must be of type 'double'.");
	
    int n = LENGTH(ss);
    int nvert = LENGTH(ssubvox);
    int k = LENGTH(smu);
	if (n <= 0)
		error("The length of 's' must be positive.");
	if (nvert <= 0)
		error("The number of voxels must be positive.");
    if (k <= 0)
		error("The number of components must be positive.");
	if (k != LENGTH(ssigma))
		error("Length of 'mu' and 'sigma' do not match.");
	int ldZ = LENGTH(sZ) / k;
    if (ldZ <= 0)
		error("The leading dimension of 'Z' must be positive.");
	int nsub = nvert / n;
    if (nvert != ldZ - 1)
		error("Number of indices and observations do not match.");
	
	
   	double *s = REAL(ss);
    int *subvox = INTEGER(ssubvox);
    int *Z = INTEGER(sZ);
    double *mu = REAL(smu);
    double *sigma = REAL(ssigma);
    SEXP val = PROTECT(allocVector(REALSXP, nvert));
    double *Y = REAL(val);
    int i, j, b;
    double *a = (double *) R_alloc(nsub, sizeof(double));
    double *V = (double *) R_alloc(nsub, sizeof(double));
    int *cls = (int *) R_alloc(nsub, sizeof(int));
	
    
    GetRNGstate();
	
    for (i = 0; i < n; i++) {
		double sumMu = 0.0;
		double sumSigma2 = 0.0;
		double sumV = 0.0;
		for (j = 0; j < nsub; j++) {
			/* use subvox entry to pick row of Z */
			/* then construct cls from 0s and 1s in Z row. */
			int m = subvox[i + j * n] - 1;
			for (b = 0; b < k; b++)
				if (Z[m + b * ldZ] == 1.0) {
					cls[j] = b;
					break;
				}
		}
		for (j = 0; j < nsub; j++) {
			double s = sigma[cls[j]];
			double s2 = s * s;
			a[j] = s2;
			sumSigma2 += s2;
			sumMu += mu[cls[j]];
		}
		for (j = 0; j < nsub; j++) {
			a[j] = a[j] / sumSigma2; 
			V[j] = norm_rand() * sigma[cls[j]];
			sumV += V[j];
		}
		for (j = 0; j < nsub; j++) {
			int m = subvox[i + j * n] - 1;
			double Ystar = V[j] - a[j] * sumV;
			Y[m] = Ystar + mu[cls[j]] + a[j] * (s[i] - sumMu);
		}
    }
	
    PutRNGstate();
	
    UNPROTECT(1);
    return val;
}
