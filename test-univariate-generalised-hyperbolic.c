#include "ambit-stochastics.h"

int main(int argc, char *argv[]) {

    int err = 0;

    ptrdiff_t n;
    double lambda;
    double alpha;
    double beta;
    double mu;
    double delta;

    int     output    = false;
    int     n_threads = 0;
    double *x         = NULL;

    gsl_rng_env_setup();
    const gsl_rng_type  *T   = gsl_rng_default;
    gsl_rng            **rng = NULL;

    char usage[] = "Usage: %s n lambda alpha beta mu delta\n";

    if (argc != 7) {
        printf(usage, argv[0]);
        err = -1;
        goto cleanup;
    }

    int ctr = 0;
    n       = atol(argv[++ctr]);
    output  = n > 0;
    n       = labs(n);
    lambda  = atof(argv[++ctr]);
    alpha   = atof(argv[++ctr]);
    beta    = atof(argv[++ctr]);
    mu      = atof(argv[++ctr]);
    delta   = atof(argv[++ctr]);

    x         = calloc(n, sizeof(double));
    n_threads = omp_get_max_threads();
    rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
    if (!x || !rng) {err = -1; goto cleanup;}
    for (int i = 0; i < n_threads; i++) {
        rng[i] = gsl_rng_alloc(T);
        gsl_rng_set(rng[i], i + 1);
    }
    
    err = generate_univariate_generalised_hyperbolic(
        n_threads, rng, lambda, alpha, beta, mu, delta, n, x);
    if (err) {
        printf("Error in generate_univariate_generalised_hyperbolic.\n");
        goto cleanup;
    }

    if (output) 
        for (ptrdiff_t j = 0; j < n; j++) 
            printf(" %f\n", x[j]);
    
  cleanup:
    free(x);
    return err;
}
