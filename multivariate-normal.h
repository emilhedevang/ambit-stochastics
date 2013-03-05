#include "ambit-stochastics.h"

int generate_multivariate_normal(int n_threads, gsl_rng **rng,
                                 int dim, double *mean, double *covariance,
                                 ptrdiff_t n, double *x);

