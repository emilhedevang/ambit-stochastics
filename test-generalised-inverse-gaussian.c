#include "ambit-stochastics.h"

int main (int argc, char *argv[]) {
  
  int    err    = 0;
  int    n      = 0;
  double lambda = 0.0;
  double chi    = 0.0;
  double psi    = 0.0;

  double *x         = NULL;
  int     n_threads = 0;

  gsl_rng_env_setup();
  const gsl_rng_type  *T   = gsl_rng_default;
  gsl_rng            **rng = NULL;

  if (argc != 5) {
    err = -1;
    printf("Usage: %s n lambda chi psi\n", argv[0]);
    goto cleanup;
  }

  n      = atoi(argv[1]);
  lambda = atof(argv[2]);
  chi    = atof(argv[3]);
  psi    = atof(argv[4]);

  x         = (double *) malloc(n * sizeof(double));
  n_threads = omp_get_max_threads();
  rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
  for (int i = 0; i < n_threads; i++) {
    rng[i] = gsl_rng_alloc(T);
    gsl_rng_set(rng[i], i + 1);
  }
  
  err = generate_generalised_inverse_gaussian(n_threads, rng, lambda, chi, psi, n, x);

  if (err) {
    printf("Error in generate_generalised_inverse_gaussian.\n");
    goto cleanup;
  }
  
  for (int i = 0; i < n; i++) 
    printf(" %f\n", x[i]);
  
 cleanup:
  for (int i = 0; i < n_threads; i++)
    gsl_rng_free(rng[i]);
  free(rng);
  free(x);
  
  return 0;
}
