#include "ambit-stochastics.h"

int main (int argc, char *argv[]) {

  int    err   = 0;
  int    dim   = 3;
  double mu[3] = {-2, 1, -5};
  double a[9]  = {2, 1, 0,
                  1, 3, 0,
                  0, 0, 0};
  double lambda[3];
  double u[9];
  int     n         = 0;
  double *x         = NULL;
  int     n_threads = 0;

  gsl_rng_env_setup();
  const gsl_rng_type  *T   = gsl_rng_default;
  gsl_rng            **rng = NULL;

  if (argc != 2) {
    err = -1;
    printf("Usage: %s n\n", argv[0]);
    goto cleanup;
  }
  
  err = spectral_decomposition(dim, a, lambda, u);
  if (err) {
    printf("Error in spectral_composition.\n");
    goto cleanup;
  }
  
  printf("lambda\n");
  for (int i = 0; i < dim; i++)
    printf(" %f", lambda[i]);
  printf("\n");

  printf("u\n");
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++)
      printf(" %f", u[i * dim + j]);
    printf("\n");
  }

  n         = atoi(argv[1]);
  x         = calloc(n * dim, sizeof(double));
  n_threads = omp_get_max_threads();
  rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
  for (int i = 0; i < n_threads; i++) {
    rng[i] = gsl_rng_alloc(T);
    gsl_rng_set(rng[i], i + 1);
  }
  
  err = generate_multivariate_normal(n_threads, rng, dim, mu, a, n, x);

  if (err) {
    printf("Error in generate_multivariate_normal.\n");
    goto cleanup;
  }
  
  printf("x\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++)
      printf(" %f", x[i * dim + j]);
    printf("\n");
  }
  
 cleanup:
  for (int i = 0; i < n_threads; i++)
    gsl_rng_free(rng[i]);
  free(rng);
  free(x);
  
  return err;
}
