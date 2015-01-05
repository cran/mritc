#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP getDenSubBias(SEXP sy, SEXP sbias, SEXP smu, SEXP ssigma)
{
    if (TYPEOF(sy) != REALSXP)
		error("'y' must be of type 'double'.");
    if (TYPEOF(sbias) != REALSXP)
		error("'bias' must be of type 'double'.");
   	if (TYPEOF(smu) != REALSXP)
		error("'mu' must be of type 'double'.");
    if (TYPEOF(ssigma) != REALSXP)
		error("'sigma' must be of type 'double'.");
	
	int nvert = LENGTH(sy);
    int k = LENGTH(smu);
	if (nvert <= 0)
		error("The number of voxels must be positive.");
    if (k <= 0)
		error("The number of components must be positive.");
	if (k != LENGTH(ssigma))
		error("Length of 'mu' and 'sigma' do not match.");
	
    SEXP val = allocMatrix(REALSXP, nvert, k);
    double *bias, *mu, *sigma, *y;
    int i, j;
    y = REAL(sy);
	bias = REAL(sbias);
    mu = REAL(smu);
    sigma = REAL(ssigma);
	
    for (j = 0; j < k; j++) {
		double *den = REAL(val) + j * nvert;
		double mu_j = mu[j];
		double sigma_j = sigma[j];
		double nc_j = M_1_SQRT_2PI / sigma_j;
#pragma omp parallel for firstprivate(mu_j, sigma_j, nc_j)
		for (i = 0; i < nvert; i++) {
			double res = (y[i] - bias[i] - mu_j) / sigma_j;
			den[i] = exp(-0.5 * res * res) * nc_j;
		}
    }
    return val;
}

