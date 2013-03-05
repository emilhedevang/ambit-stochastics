#include "discrete-convolutions.h"

#define WRAP(i, n) ((i) >= (n) ? (i) - (n) : (i))

/* void multidimensional_index(ptrdiff_t ind, int rank, int *sub, int *n) { */
/*     long p = ind; */
/*     for (int i = rank - 1; i >= 0; i--) { */
/*         ldiv_t qr = ldiv(p, n[i]); */
/*         if (qr.rem < 0) { */
/*             qr.quot--; */
/*             qr.rem += n[i]; */
/*         } */
/*         sub[i] = qr.rem; */
/*         p = qr.quot; */
/*     } */
/* } */

ptrdiff_t linear_index(int rank, int *dim, int *dimembed, int *sub) {
    ptrdiff_t index = 0;
    for (int i = 0; i < rank; i++) 
        index = index * dimembed[i] + WRAP(sub[i], dim[i]);
    return index;
}

void increment_subscript(int rank, int *dim, int *sub) {
    for (int i = rank - 1; i >= 0; i--)
        if (++sub[i] < dim[i])
            return;
        else 
            sub[i] = 0;
}

void embed_array_double(int rank, 
                        int *xdim, int *xdimembed, double *x, 
                        int *ydim, int *ydimembed, double *y, 
                        int *offset) {
    int *xsub = NULL;
    int *ysub = NULL;
    ptrdiff_t xindex = 0;
    ptrdiff_t yindex = 0;
    ptrdiff_t xn     = 1;
    for (int i = 0; i < rank; i++) 
        xn *= xdim[i];

#pragma omp parallel private(xsub, ysub, xindex, yindex)
    {
        xsub = (int *)malloc(rank * sizeof(int));
        ysub = (int *)malloc(rank * sizeof(int));
        
#pragma omp for
        for (int i0 = 0; i0 < xdim[0]; i0++) {
            /* Setup the first subscript. */
            xsub[0] = i0;
            for (int i = 1; i < rank; i++)
                xsub[i] = 0;
            
            /* Loop over the rest of the subscripts. */
            for (ptrdiff_t i = 0; i < xn / xdim[0]; i++) {
                for (int j = 0; j < rank; j++)
                    ysub[j] = xsub[j] + offset[j];
                xindex = linear_index(rank, xdim, xdimembed, xsub);
                yindex = linear_index(rank, ydim, ydimembed, ysub);
                y[yindex] = x[xindex];
                increment_subscript(rank, xdim, xsub);
            }
        }
    }
}

/** 
 * Calculates the circular convolution of `k` and `x`, overwriting `x`
 * with the result.
 *
 * The array `x` is assumed to be padded along the last dimension with
 * 2 extra elements if `xdim[rank - 1]` is even and 1 extra element if
 * `xdim[rank - 1]` is odd. For example, 
 * 
 *     [ * * * * * * 0 0 ]
 *     [ * * * * * * 0 0 ]
 *     [ * * * * * * 0 0 ]
 *
 * and
 *
 *     [ * * * * * 0 ]
 *     [ * * * * * 0 ]
 *     [ * * * * * 0 ]
 * 
 */
int ambit_circular_convolution_inplace(int rank, 
                                       int *kdim, double *k, 
                                       int *xdim, double *x) {
    int err = AMBIT_NO_ERROR;

    int n_threads = omp_get_max_threads();
    fftw_plan_with_nthreads(n_threads);
    
    /*
     * Decleare and variables whose memory is to be freed during
     * cleanup.
     */
    
    double *y = NULL;
    fftw_plan x_plan = NULL;
    fftw_plan y_plan = NULL;
    fftw_plan X_plan = NULL;

    /*
     * Calculate the implicit embedding dimensions of `x`.
     */
    int *xdimembed = (int *)malloc(rank * sizeof(int));
    for (int i = 0; i < rank - 1; i++)
        xdimembed[i] = xdim[i];
    xdimembed[rank - 1] = 2 * (xdim[rank - 1] / 2 + 1);
                     
    /* 
     * Embed `k` in the larger array `y` initialized with zeroes.
     * The array `y` has the same dimensions as the array `x`.
     */
    ptrdiff_t xn = 1;
    ptrdiff_t xnembed = 1;
    for (int i = 0; i < rank; i++) {
        xn *= xdim[i];
        xnembed *= xdimembed[i];
    }
    
    y = (double *)fftw_malloc(xnembed * sizeof(double));
    if (!y) {
        err = AMBIT_ALLOCATION_ERROR; 
        goto cleanup;
    }
    
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < xnembed; i++) 
        y[i] = 0.0;

    int *offset = (int *)malloc(rank * sizeof(int));
    for (int i = 0; i < rank; i++)
        offset[i] = xdim[i] - kdim[i] + 1;
    
    embed_array_double(rank, kdim, kdim, k, xdim, xdimembed, y, offset);

    /*
     * Plan and calculate DTFs
     */
    fftw_complex *X = (fftw_complex *)x;
    fftw_complex *Y = (fftw_complex *)y;
    x_plan = fftw_plan_dft_r2c(rank, xdim, x, X, FFTW_ESTIMATE);
    y_plan = fftw_plan_dft_r2c(rank, xdim, y, Y, FFTW_ESTIMATE);
    if (!x_plan || !y_plan) {
        err = AMBIT_PLAN_ERROR;
        goto cleanup;
    }
    fftw_execute(x_plan);
    fftw_execute(y_plan);

    /*
     * Multiply DFTs
     */
    double scale = 1 / (double)xn;
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < xnembed / 2; i++)
        X[i] *= scale * Y[i];

    /*
     * Plan and calculate IDFT
     */
    X_plan = fftw_plan_dft_c2r(rank, xdim, X, x, FFTW_ESTIMATE);
    if (!X_plan) {
        err = AMBIT_PLAN_ERROR;
        goto cleanup;
    }
    fftw_execute(X_plan);

    /*
     * Cleanup and return
     */
cleanup:
    fftw_destroy_plan(y_plan);
    fftw_destroy_plan(x_plan);
    fftw_destroy_plan(X_plan);
    fftw_free(y);

    return err;    
}

/**
 * Calculates the circular matrix convolution of `k` and `x` in place,
 * overwriting `x` with the result.
 *
 * Some explanation is required. Let `K` denote a matrix with `nrows`
 * rows and `ncols`. Each entry `K(i,j)` is a scalar array of rank
 * `rank` and with dimensions given by `kdim`. Let `X` denote a vector
 * of length `ncols`. Each entry `X(j)` is a scalar array of rank
 * `rank` and dimensions given by `xdim`. By the circular matrix
 * convolution of `K` and `X` we understand a vector `Y` of length
 * `nrows` whose `i`'th entry `Y(i)` is a scalar array of rank `rank`
 * and dimensions given by `xdim`, and `Y(i)` is equal to the sum of
 * the scalar circular convolutions of `K(i,j)` and `X(j)` for `j = 0,
 * ..., kncols - 1`.
 *
 * `X` is represented in memory by the arrays starting at the memory
 * locations `x[0]`, ..., `x[ncols - 1]`. Each array `x[i]` is assumed
 * to be padded as described for
 * `ambit_circular_convolution_inplace`. The array behind `K(i,j)`
 * starts at `k[matindex[i * ncols + j]]`. Thus, `matindex` is an
 * array of length `nrows * ncols`.
 */
int ambit_circular_matrix_convolution_inplace(int rank,
                                              int *kdim, double **k,
                                              int *xdim, double **x, 
                                              int nrows, int ncols, int *matindex) {
    int err = AMBIT_NO_ERROR;

    /*
     * Inplace DFT of x
     */
    fftw_plan x_plan = NULL;
    fftw_complex **X = (fftw_complex **)malloc(ncols * sizeof(fftw_complex *));
    for (int i = 0; i < ncols; i++) {
        X[i] = (fftw_complex)x[i];
        x_plan = fftw_plan_dft_r2c(rank, xdim, x, X, FFTW_ESTIMATE);
        if (!x_plan) {
            err = AMBIT_PLAN_ERROR;
            goto cleanup;
        }
        fftw_execute(x_plan);
        fftw_destroy_plan(x_plan);
    }
    
    /*
     * Allocate temporary storage
     */
    ptrdiff_t xnembed = 1;
    for (int i = 0; i < rank - 1; i++) 
        xnembed *= xdim[i];
    xnembed *= 2 * (xdim[rank - 1] / 2 + 1);
            
    double       *t  = NULL;
    fftw_complex *T  = NULL;
    t = (double *)fftw_malloc(xnembed * sizeof(double));
    T = (fftw_complex *)t;
    
    fftw_complex **Y = (fftw_complex **)malloc(nrows * sizeof(fftw_complex *));
    for(int i = 0; i < nrows; i++)
        Y[i] = (fftw_complex *)fftw_malloc(xnembed / 2 * sizeof(fftw_complex));
        

    /*
     * Cleanup and return
     */
cleanup:
    fftw_destroy_plan(x_plan);
    


    return err;
}
