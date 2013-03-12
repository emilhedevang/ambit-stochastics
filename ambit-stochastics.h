#ifndef AMBIT_STOCHASTICS_H
#define AMBIT_STOCHASTICS_H

/* Standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

/* GNU Scientific Library */
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* OpenMP */
#include <omp.h>

/* LAPACKE and BLAS */
#include <lapacke.h>
#include <atlas/cblas.h>

/* HDF 5 */
#include "hdf5.h"
#include "hdf5_hl.h"


/********************* 
 * Ambit Stochastics *
 *********************/

/*
 * generalised-inverse-gaussian
 */

int generate_generalised_inverse_gaussian(int n_threads, gsl_rng **rng,
                                          double lambda, double chi, double psi, 
                                          ptrdiff_t n, double *x);

/*
 * multivariate-generalised-hyperbolic
 */

int generate_multivariate_generalised_hyperbolic(
    int n_threads, gsl_rng **rng,
    int dim,
    double lambda, double alpha, double *beta, double *mu, double delta, double *Delta,
    ptrdiff_t n, double *x);

int generate_univariate_generalised_hyperbolic(
    int n_threads, gsl_rng **rng,
    double lambda, double alpha, double beta, double mu, double delta,
    ptrdiff_t n, double *x);

/*
 * multivariate-normal
 */

int generate_multivariate_normal(int n_threads, gsl_rng **rng,
                                 int dim,
                                 double *mean, double *covariance,
                                 ptrdiff_t n, double *x);

int generate_univariate_normal(int n_threads, gsl_rng **rng,
                               double mean, double variance,
                               ptrdiff_t n, double *x);
/* 
 * utilities
 */

int spectral_decomposition(int dim, double *a, double *lambda, double *u);


/*
 * trawl-process
 */

typedef struct {double mean, variance;} 
    normal_generator_params;

typedef struct {double lambda, alpha, beta, mu, delta;} 
    univariate_generalised_hyperbolic_generator_params;

int constant_generator(int n_threads, gsl_rng **rng, 
                       double volume, void *p, 
                       int64_t n, double *x);

int normal_generator(int n_threads, gsl_rng **rng, 
                     double volume, void *p, 
                     int64_t n, double *x);

int univariate_generalised_hyperbolic_generator(int n_threads, gsl_rng **rng, 
                                                double volume, void *p, 
                                                int64_t n, double *x);

int trawl_process(
    int n_threads, gsl_rng **rng,
    int (*generator)(int, gsl_rng **, double volume, void *, int64_t, double *),
    void *generator_params,
    int64_t n_slices, double *slice_volume,
    int64_t n, double *x);


#endif

