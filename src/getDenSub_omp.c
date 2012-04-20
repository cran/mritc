#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP getDenSub(SEXP sy, SEXP smu, SEXP ssigma)
{
    if (TYPEOF(sy) != REALSXP)
		error("'y' must be of type 'double'.");
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
    double *mu, *sigma, *y;
    int i, j;
    y = REAL(sy);
    mu = REAL(smu);
    sigma = REAL(ssigma);
	
    for (j = 0; j < k; j++) {
		double *den = REAL(val) + j * nvert;
		double mu_j = mu[j];
		double sigma_j = sigma[j];
		double nc_j = M_1_SQRT_2PI / sigma_j;
#pragma omp parallel for firstprivate(mu_j, sigma_j, nc_j)
		for (i = 0; i < nvert; i++) {
			double res = (y[i] - mu_j) / sigma_j;
			den[i] = exp(-0.5 * res * res) * nc_j;
		}
    }
    return val;
}

#define TABSIZ 20000
static double dnrm[TABSIZ * 3];

SEXP getDenSubTable(SEXP sy, SEXP snvert, SEXP smu, SEXP ssigma, SEXP sk)
{
    int nvert = asInteger(snvert);
    int k = asInteger(sk);
    SEXP val = allocMatrix(REALSXP, nvert, k);
    double *mu, *sigma, *y;
    double ymin, ymax;
    int i, j;
	
    /**** error checks */
    y = REAL(sy);
    mu = REAL(smu);
    sigma = REAL(ssigma);
	
    ymin = ymax = y[0];
    for (i = 1; i < nvert; i++)
		if (y[i] < ymin) ymin = y[i];
		else if (y[i] > ymax) ymax = y[i];
	
    for (j = 0; j < k; j++) {
		double mu_j = mu[j];
		double sigma_j = sigma[j];
		double nc_j = M_1_SQRT_2PI / sigma_j;
		for (i = 0; i < TABSIZ; i++) {
			double x = ymin + (i * (ymax - ymin)) / (TABSIZ - 1);
			double res = (x - mu_j) / sigma_j;
			dnrm[i + j * TABSIZ] = exp(-0.5 * res * res) * nc_j;
		}
    }
	
	
    for (j = 0; j < k; j++) {
		double *den = REAL(val) + j * nvert;
		/*double h = (TABSIZ - 1) / (ymax - ymin);*/
		for (i = 0; i < nvert; i++) {
			int m = (TABSIZ - 1) * (y[i] - ymin) / (ymax - ymin);
			/*	    if (m == TABSIZ - 1)
					den[i] = dnrm[m + j * TABSIZ];
					else {
					double d1 = dnrm[m + j * TABSIZ];
					double d2 = dnrm[m + 1 + j * TABSIZ];
					double p = (y[i] - ymin - (m * (ymax - ymin)) / (TABSIZ - 1)) / h;
					den[i] = p * d2 + (1 - p) * d1;
					}*/
			den[i] = dnrm[m + j * TABSIZ];
		}
    }
    return val;
}

