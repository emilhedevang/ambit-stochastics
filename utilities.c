#include "utilities.h"

int spectral_decomposition(int dim, double *a, double *lambda, double *u) {
  int err              = 0;
  lapack_int  info     = 0;
  lapack_int  ev_found = 0;
  double     *a_copy   = NULL;
  lapack_int *isuppz   = NULL;

  if (dim < 1 || !a || !lambda || !u) {
    err = -1;
    goto cleanup;
  }

  a_copy = calloc(dim * dim, sizeof(double));
  isuppz = calloc(2 * dim  , sizeof(lapack_int));
  
  if (!a_copy || !isuppz) {
    err = -1;
    goto cleanup;
  }
  
  cblas_dcopy(dim * dim, a, 1, a_copy, 1);
  info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'L', 
                        dim, a_copy, dim, 
                        0.0, 0.0, 0, 0, 0.0, 
                        &ev_found, lambda, u, dim, isuppz);

  if (info || ev_found < dim) {
    err = -1;
    goto cleanup;
  }

 cleanup:
  free(a_copy);
  free(isuppz);
  return err;
}
