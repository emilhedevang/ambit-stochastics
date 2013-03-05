#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include <omp.h>

#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include <atlas/cblas.h>

enum ambit_error_code {
    AMBIT_NO_ERROR         =  0,
    AMBIT_INVALID_INPUT    = -1,
    AMBIT_ALLOCATION_ERROR = -2,
    AMBIT_PLAN_ERROR       = -3
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
