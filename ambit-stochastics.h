#ifndef AMBIT_STOCHASTICS_H
#define AMBIT_STOCHASTICS_H

/* Standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <complex.h>
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

/* FFTW 3 */
#include <fftw3.h>


/********************* 
 * Ambit Stochastics *
 *********************/

/* 
 * discrete-convolution
 */

enum ambit_error_code {
    AMBIT_NO_ERROR         =  0,
    AMBIT_INVALID_INPUT    = -1,
    AMBIT_ALLOCATION_ERROR = -2,
    AMBIT_FFTW_PLAN_ERROR  = -3
};



ptrdiff_t linear_index(int rank, int *dim, int *dimembed, int *sub);

void embed_array_double(int rank, 
                        int *xdim, int *xdimembed, double *x, 
                        int *ydim, int *ydimembed, double *y, 
                        int *offset);

void increment_subscript(int rank, int *dim, int *sub);

int ambit_circular_convolution_inplace(int rank, 
                                       int *kdim, double *k, 
                                       int *xdim, double *x);

int ambit_symmetric_odd_isotropic_circular_convolution_inplace(
    int n_threads,
    int n_k, double *k_abscissa, double *k_ordinate,
    int dim, double delta_x, double *x1, double *x2, double *x3);


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

void print_array_double_3(int *dims, double *x);

double linear_interpolation(double x, double a, double b, int n, double *tbl, double y0);

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

