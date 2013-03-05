#include "ambit-stochastics.h"

int generate_generalised_inverse_gaussian(int n_threads, gsl_rng **rng,
                                          double lambda, double chi, double psi, 
                                          ptrdiff_t n, double *x);
