#include "multivariate-generalised-hyperbolic.h"

int generate_multivariate_generalised_hyperbolic(
    int n_threads, gsl_rng **rng,
    int dim,
    double lambda, double alpha, double *beta, double *mu, double delta, double *Delta,
    ptrdiff_t n, double *x) {
    
    int     err        = 0;
    double *w          = NULL;
    double *Delta_beta = NULL;
    
    w          = calloc(n, sizeof(double));
    Delta_beta = calloc(dim, sizeof(double));
    if (!w || !Delta_beta) { err = -1; goto cleanup;}
    
    /*
     * Generate GIG 
     */
    double chi = delta * delta;
    double psi;
    cblas_dsymv(CblasColMajor, CblasLower, dim, 1.0, Delta, dim, beta, 1, 0.0, Delta_beta, 1);
    psi = alpha * alpha - cblas_ddot(dim, beta, 1, Delta_beta, 1);
    err = generate_generalised_inverse_gaussian(n_threads, rng, lambda, chi, psi, n, w);
    if (err) goto cleanup;
    
    /*
     * Generate normal
     */
    err = generate_multivariate_normal(n_threads, rng, dim, NULL, Delta, n, x);
    if (err) goto cleanup;
    
    /* 
     * Combine to GH 
     */
#pragma omp parallel for num_threads(n_threads)
    for (ptrdiff_t j = 0; j < n; j++) {
        cblas_dscal(dim, sqrt(w[j]), x + j * dim, 1);
        cblas_daxpy(dim, w[j], Delta_beta, 1, x + j * dim, 1);
        cblas_daxpy(dim, 1.0, mu, 1, x + j * dim, 1);
    }
    
  cleanup:
    free(w);
    free(Delta_beta);
    return err;
}


