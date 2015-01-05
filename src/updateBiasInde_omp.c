#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Destructively updates the bias field in sbias */
/* openMP is used*/
SEXP updateBiasIndeX(SEXP sy, SEXP sbias, SEXP sneighbors, SEXP snneigh, SEXP snneighbors, SEXP sZ, SEXP smu, SEXP ssigma, SEXP sgamma, SEXP sweights)
{
	if (TYPEOF(sy) != REALSXP)
		error("'y' must be of type 'double'.");
	if (TYPEOF(sbias) != REALSXP)
		error("'bias' must be of type 'double'.");
	if (TYPEOF(sneighbors) != INTSXP)
		error("'neighbors' must be of type 'integer'.");
	if (TYPEOF(snneigh) != INTSXP)
		error("'nneigh' must be of type 'integer'.");
	if (TYPEOF(snneighbors) != REALSXP)
		error("'nneighbors' must be of type 'double'.");
    if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'.");
 	if (TYPEOF(smu) != REALSXP)
		error("'mu' must be of type 'double'.");
 	if (TYPEOF(ssigma) != REALSXP)
		error("'sigma' must be of type 'double'.");
 	if (TYPEOF(sgamma) != REALSXP)
		error("'gamma' must be of type 'double'.");
 	if (TYPEOF(sweights) != REALSXP)
		error("'weights' must be of type 'double'.");

    double *y = REAL(sy);
    double *bias = REAL(sbias);
	int *neighbors = INTEGER(sneighbors);
    int nneigh = asInteger(snneigh);
    if (nneigh <= 0)
		error("The number of neighbors must be positive.");
	double *nneighbors = REAL(snneighbors);
	int k = LENGTH(smu);
	if (k <= 0)
		error("The number of components must be positive.");
	int *Z = INTEGER(sZ);
    int ldZ = LENGTH(sZ) / k;
    if (ldZ <= 0)
		error("The leading dimension of 'Z' must be positive.");
	int ldN = LENGTH(sneighbors) / nneigh;
	if (ldN <= 0)
		error("The leading dimension of 'neighbors' must be positive.");
	if (ldZ - 1 != ldN || ldN != LENGTH(snneighbors))
		error("The leading dimension of 'Z' and 'neighbors' and the length of 'nneighbors' do not match.");
	double *mu = REAL(smu);
	double *sigma = REAL(ssigma);
	double gamma = asReal(sgamma);
	double *weights = REAL(sweights);

	int n = LENGTH(snneighbors);
	int i;

#pragma omp parallel for firstprivate(k, ldN, ldZ, nneigh, neighbors, \
				      nneighbors, Z, weights)
	for (i = 0; i < n; i++) {
		int m, j, ind=0;
		double sumB = 0.0;
		double meanN, preN; 
		
		/* compute the sum of bias of neighbors */
		for (m = 0; m < nneigh; m++) {
			int mm = neighbors[i + m * ldN] - 1;
			sumB = sumB + bias[mm]*weights[m];
		}

		/* get the tissue type of voxel i */
		for (j = 0; j < k; j++) {
			if (Z[i + j * ldZ] == 1) {
				ind = j;
			}
		}
	
		double ni = nneighbors[i];
		double tau = 1 / (sigma[ind] * sigma[ind]);
		preN = ni * gamma + tau;
		meanN = gamma * sumB + tau * (y[i] - mu[ind]);
		meanN = meanN / preN;
		bias[i] =  meanN;
	}

	return sbias;
}

/* Non-destructive version */
SEXP updateBiasInde(SEXP sy, SEXP sbias, SEXP sneighbors, SEXP snneigh, SEXP snneighbors, SEXP sZ, SEXP smu, SEXP ssigma, SEXP sgamma, SEXP sweights)
{
	SEXP val;
    PROTECT(sbias = duplicate(sbias));
    val = updateBiasIndeX(sy, sbias, sneighbors, snneigh, snneighbors, sZ,
					  smu, ssigma, sgamma, sweights);
    UNPROTECT(1);
    return val;
}
