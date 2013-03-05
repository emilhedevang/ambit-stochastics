#include "multivariate-normal.h"

/** Generate multivariate normal variates.
 *
 * n_threads: number of threads to use in OpenMP
 * rng: random number generators, array of length n_threads
 * dim: number of dimensions
 * mu: mean value, array of length dim
 * sigma: covariance matrix, vector of length dim*dim
 * n: number of variates to generate
 * x: storage of the result, array of length n*dim
 */
int generate_multivariate_normal(int n_threads, gsl_rng **rng,
                                 int dim, double *mean, double *covariance,
                                 ptrdiff_t n, double *x) {
  int err        = 0;
  double *U      = NULL;
  double *lambda = NULL;
  double *tmp    = NULL;

  if (n_threads < 1 || dim < 1 || n < 1 || !rng || !covariance || !x) {
    err = -1;
    goto cleanup;
  }

  U      = calloc(dim * dim, sizeof(double));
  lambda = calloc(dim      , sizeof(double));
  tmp    = calloc(dim      , sizeof(double));
  
  if (!U || !lambda || !tmp) {
    err = -1;
    goto cleanup;
  }

  err = spectral_decomposition(dim, covariance, lambda, U);
  if (err) goto cleanup;
  
  for (int i = 0; i < dim; i++)
    if (lambda[i] < 0.0) {
      err = -1;
      goto cleanup;
    }
  
  double   sq        = 0.0;
  int      thread_id = 0;
  gsl_rng *rng0      = NULL;

  for (int i = 0; i < dim; i++) {
    sq = sqrt(lambda[i]);
    
#pragma omp parallel private(thread_id, rng0) num_threads(n_threads)
    {
      thread_id = omp_get_thread_num();
      rng0 = rng[thread_id];
      
#pragma omp for
      for (ptrdiff_t j = 0; j < n; j++)
        x[i + dim * j] = sq * gsl_ran_ugaussian(rng0);
    }
  }
  
  /* Should be expressed in terms of dgemm when it supports 64 bit matrix sizes */
  /* Can be accellerated using OpenMP */
  for (ptrdiff_t j = 0; j < n; j++) {
    cblas_dcopy(dim, x + j * dim, 1, tmp, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, 
                dim, dim, 1.0, U, dim, tmp, 1, 0.0, x + j * dim, 1);
  }

  /* Handle non-zero mean */
  if (mean)
      for (ptrdiff_t j = 0; j < n; j++) {
          cblas_daxpy(dim, 1.0, mean, 1, x + j * dim, 1);
      }
  
 cleanup:
  free(U);
  free(lambda);
  return err;
}

