#define _XOPEN_SOURCE 700

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


//
// Dense arrays
//

#define AMBIT_ALIGNMENT 32

typedef double R;
typedef fftw_complex C;

struct ambit_dense_array {
    // Primary fields
    unsigned int rank;
    size_t *dim;
    size_t *dimembed;
    R      *data;
    // Derived fields
    size_t n;            // product of dim
    size_t nembed;       // product of dimembed
    size_t *strideembed; // derived from dimembed
};

size_t
wrap_subscript(size_t sub, size_t dim);

void
increment_subscript(size_t *sub, int rank, size_t *dim);

struct ambit_dense_array *
ambit_dense_array_malloc (int rank, size_t *dim, size_t *dimembed);

struct ambit_dense_array *
ambit_dense_array_malloc_1d(size_t dim0);

struct ambit_dense_array *
ambit_dense_array_malloc_2d(size_t dim0, size_t dim1);

struct ambit_dense_array *
ambit_dense_array_malloc_3d(size_t dim0, size_t dim1, size_t dim2);

void 
ambit_dense_array_free (struct ambit_dense_array *array);

bool 
ambit_dense_array_valid (struct ambit_dense_array *a);

bool
ambit_dense_array_ranks_eq (struct ambit_dense_array *a, struct ambit_dense_array *b);

bool
ambit_dense_array_dims_eq (struct ambit_dense_array *a, struct ambit_dense_array *b);

bool
ambit_dense_array_dims_leq (struct ambit_dense_array *a, struct ambit_dense_array *b);

bool 
ambit_dense_array_fftw_r2c_embedded_tightly(struct ambit_dense_array *a);

bool 
ambit_dense_array_fftw_r2c_embedded(struct ambit_dense_array *a);

void
ambit_dense_array_fprintf (FILE *stream, struct ambit_dense_array *array, bool also_data);

size_t
ambit_dense_array_linear_index(struct ambit_dense_array *array, size_t *sub);

size_t
ambit_dense_array_linear_index_1d(struct ambit_dense_array *array, size_t sub0);

size_t
ambit_dense_array_linear_index_2d(struct ambit_dense_array *array, size_t sub0, size_t sub1);

size_t
ambit_dense_array_linear_index_3d(struct ambit_dense_array *array, 
                                  size_t sub0, size_t sub1, size_t sub2);

void 
ambit_dense_array_set_all(struct ambit_dense_array *a, R c);

void
ambit_dense_array_set(struct ambit_dense_array *array, size_t *sub, R value);

int
ambit_dense_array_circular_embed (struct ambit_dense_array *a,
                                  struct ambit_dense_array *b,
                                  size_t *offset);
void
ambit_dense_array_parallel_fill(struct ambit_dense_array *a,
                                R (*func)(const void *, const void *),
                                const void  *p,  size_t p_size,
                                const void **pp, size_t pp_size);

//
// Discrete Fourier Transforms
//

int
ambit_dft_r2c_inplace(struct ambit_dense_array *a, int sign);

int 
ambit_circular_convolution_inplace(struct ambit_dense_array *k, struct ambit_dense_array *x);





/* 
 * discrete-convolution
 */

enum ambit_error_code {
    AMBIT_NO_ERROR         =  0,
    AMBIT_INVALID_INPUT    = -1,
    AMBIT_ALLOCATION_ERROR = -2,
    AMBIT_FFTW_PLAN_ERROR  = -3
};



/* ptrdiff_t linear_index(int rank, int *dim, int *dimembed, int *sub); */

/* void embed_array_double(int rank,  */
/*                         int *xdim, int *xdimembed, double *x,  */
/*                         int *ydim, int *ydimembed, double *y,  */
/*                         int *offset); */

//void increment_subscript(int rank, int *dim, int *sub);

/* int ambit_circular_convolution_inplace(int rank,  */
/*                                        int *kdim, double *k,  */
/*                                        int *xdim, double *x); */

/* int ambit_symmetric_odd_isotropic_circular_convolution_inplace( */
/*     int n_threads, */
/*     int n_k, double *k_abscissa, double *k_ordinate, */
/*     int dim, double delta_x, double *x1, double *x2, double *x3); */


/*
 * generalised-inverse-gaussian
 */

int 
generate_generalised_inverse_gaussian (gsl_rng **rng, struct ambit_dense_array *a, R lambda, R chi, R psi);


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

int 
generate_univariate_normal(gsl_rng **rng, struct ambit_dense_array *a, R mean, R variance);

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

