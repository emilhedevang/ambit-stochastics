#include "ambit-stochastics.h"

int generate_multivariate_generalised_hyperbolic(
    int n_threads, gsl_rng **rng,
    int dim,
    double lambda, double alpha, double *beta, double *mu, double delta, double *Delta,
    ptrdiff_t n, double *x);
