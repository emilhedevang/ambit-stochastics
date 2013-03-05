#include "ambit-stochastics.h"

int main(int argc, char *argv[]) {

    int err = 0;
    int dim = 0;

    int n;
    double lambda;
    double alpha;
    double *beta;
    double *mu;
    double delta;
    double *Delta;

    int     n_threads = 0;
    double *x = NULL;

    gsl_rng_env_setup();
    const gsl_rng_type  *T   = gsl_rng_default;
    gsl_rng            **rng = NULL;

    char usage[] = "Usage: %s n dim lambda alpha beta_1 ... beta_dim mu_1 ... mu_dim delta Delta_{1,1} ... Delta_{dim, dim}\n";

    if (argc <= 3) {
        printf(usage, argv[0]);
        err = -1;
        goto cleanup;
    }

    int ctr = 0;
    n      = atoi(argv[++ctr]);
    dim    = atoi(argv[++ctr]);
    
    if (argc != 3 + 1 + 1 + dim + dim + 1 + dim * dim) {
        printf(usage, argv[0]);
        err = -1;
        goto cleanup;
    }

    beta  = calloc(dim, sizeof(double));
    mu    = calloc(dim, sizeof(double));
    Delta = calloc(dim * dim, sizeof(double));
    if (!beta || !mu || !Delta) {err = -1; goto cleanup;}
    
    lambda = atof(argv[++ctr]);
    alpha  = atof(argv[++ctr]);
    for (int i = 0; i < dim; i++)
        beta[i] = atof(argv[++ctr]);
    for (int i = 0; i < dim; i++)
        mu[i] = atof(argv[++ctr]);
    delta  = atof(argv[++ctr]);
    for (int i = 0; i < dim * dim; i++)
        Delta[i] = atof(argv[++ctr]);

    x = calloc(n * dim, sizeof(double));
    n_threads = omp_get_max_threads();
    rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
    if (!x || !rng) {err = -1; goto cleanup;}
    for (int i = 0; i < n_threads; i++) {
        rng[i] = gsl_rng_alloc(T);
        gsl_rng_set(rng[i], i + 1);
    }
    
    err = generate_multivariate_generalised_hyperbolic(
        n_threads, rng,
        dim, lambda, alpha, beta, mu, delta, Delta,
        n, x);
    if (err) {
        printf("Error in generate_multivariate_generalised_hyperbolic.\n");
        goto cleanup;
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < dim; i++)
            printf(" %f", x[i + j * dim]);
        printf("\n");
    }
 
  cleanup:
    free(x);
    return err;
}
