#include "ambit-stochastics.h"

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
        err = AMBIT_FFTW_PLAN_ERROR;
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
        err = AMBIT_FFTW_PLAN_ERROR;
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
 * Computes
 * 
 *     [   0  k3 -k2 ]   [ x1 ]   [  k3*x2 - k2*x3 ]
 *     [ -k3   0  k1 ] * [ x2 ] = [ -k3*x1 + k1*x3 ]
 *     [  k2 -k1   0 ]   [ x3 ]   [  k2*x1 - k1*x2 ]
 *
 * and stores the result in x1, x2, x3. It is assumed that all six
 * arrays have the same dimensions and that they are padded to allow
 * in-place FFTs.
 */
int ambit_odd_isotropic_circular_convolution_inplace(
    int n_threads,
    int dim0, int dim1, int dim2,
    double *k1, double *k2, double *k3,
    double *x1, double *x2, double *x3) {

    int err = AMBIT_NO_ERROR;
    size_t nembed = (size_t)dim0 * (size_t)dim1 * (size_t)(2 * (dim2 / 2 + 1)); 
    
    fftw_plan_with_nthreads(n_threads);
    
    /*
     * Declare variables whose memory is to be freed during cleanup.
     */

    fftw_plan k1_plan = NULL;
    fftw_plan k2_plan = NULL;
    fftw_plan k3_plan = NULL;

    fftw_plan x1_plan = NULL;
    fftw_plan x2_plan = NULL;
    fftw_plan x3_plan = NULL;

    fftw_plan y1_plan = NULL;
    fftw_plan y2_plan = NULL;
    fftw_plan y3_plan = NULL;
                     
    /* /\*  */
    /*  * Allocate temporary storage  */
    /*  *\/ */
    /* y1 = (double *)fftw_malloc(nembed * sizeof(double)); */
    /* y2 = (double *)fftw_malloc(nembed * sizeof(double)); */
    /* y3 = (double *)fftw_malloc(nembed * sizeof(double)); */
    /* if (!y1 || !y2 || !y3) { */
    /*     err = AMBIT_ALLOCATION_ERROR;  */
    /*     goto cleanup; */
    /* } */
    
    /*
     * Plan DFTs and IDFTs
     */
    fftw_complex *X1 = (fftw_complex *)x1;
    fftw_complex *X2 = (fftw_complex *)x2;
    fftw_complex *X3 = (fftw_complex *)x3;

    fftw_complex *K1 = (fftw_complex *)k1;
    fftw_complex *K2 = (fftw_complex *)k2;
    fftw_complex *K3 = (fftw_complex *)k3;

    k1_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, k1, K1, FFTW_ESTIMATE);
    k2_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, k2, K2, FFTW_ESTIMATE);
    k3_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, k3, K3, FFTW_ESTIMATE);

    x1_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, x1, X1, FFTW_ESTIMATE);
    x2_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, x2, X2, FFTW_ESTIMATE);
    x3_plan = fftw_plan_dft_r2c_3d(dim0, dim1, dim2, x3, X3, FFTW_ESTIMATE);

    y1_plan = fftw_plan_dft_c2r_3d(dim0, dim1, dim2, X1, x1, FFTW_ESTIMATE);
    y2_plan = fftw_plan_dft_c2r_3d(dim0, dim1, dim2, X2, x2, FFTW_ESTIMATE);
    y3_plan = fftw_plan_dft_c2r_3d(dim0, dim1, dim2, X3, x3, FFTW_ESTIMATE);

    if (   !k1_plan || !k2_plan || !k3_plan 
        || !x1_plan || !x2_plan || !x3_plan 
        || !y1_plan || !y2_plan || !y3_plan) {
        err = AMBIT_FFTW_PLAN_ERROR;
        goto cleanup;
    }

    /*
     * Execute DFTs
     */
    fftw_execute(k1_plan);
    fftw_execute(k2_plan);
    fftw_execute(k3_plan);

    fftw_execute(x1_plan);
    fftw_execute(x2_plan);
    fftw_execute(x3_plan);

    /*
     * Multiply DFTs
     */
    fftw_complex Y1, Y2, Y3;
    double scale = 1 / (double)(dim0 * dim1 * dim2);
#pragma omp parallel for private(Y1, Y2, Y3)
    for (ptrdiff_t i = 0; i < nembed / 2; i++) {
        Y1 = scale * (K3[i] * X2[i] - K2[i] * X3[i]);
        Y2 = scale * (K1[i] * X3[i] - K3[i] * X1[i]);
        Y3 = scale * (K2[i] * X1[i] - K1[i] * X2[i]);
        X1[i] = Y1;
        X2[i] = Y2;
        X3[i] = Y3;
    }
    
    /*
     * Execute IDFTs
     */
    fftw_execute(y1_plan);
    fftw_execute(y2_plan);
    fftw_execute(y3_plan);

    /*
     * Cleanup and return
     */
cleanup:
    fftw_destroy_plan(k1_plan);
    fftw_destroy_plan(k2_plan);
    fftw_destroy_plan(k3_plan);
    fftw_destroy_plan(x1_plan);
    fftw_destroy_plan(x2_plan);
    fftw_destroy_plan(x3_plan);
    fftw_destroy_plan(y1_plan);
    fftw_destroy_plan(y2_plan);
    fftw_destroy_plan(y3_plan);
    return err;    
}





/**
 * Computes
 * 
 *     [   0  k3 -k2 ]   [ x1 ]   [  k3*x2 - k2*x3 ]
 *     [ -k3   0  k1 ] * [ x2 ] = [ -k3*x1 + k1*x3 ]
 *     [  k2 -k1   0 ]   [ x3 ]   [  k2*x1 - k1*x2 ]
 *
 * and stores the result in x1, x2, x3. It is assumed that all six
 * arrays have the same dimensions and that they are padded to allow
 * in-place FFTs.
 *
 * Assumes that dim is even.
 */
int ambit_symmetric_odd_isotropic_circular_convolution_inplace(
    int n_threads,
    int n_k, double *k_abscissa, double *k_ordinate,
    int dim, double delta_x, double *x1, double *x2, double *x3) {
    
    int err = AMBIT_NO_ERROR;

    int last_dim  = 2 * (dim / 2 + 1);
    size_t nembed = (size_t)dim * (size_t)dim * (size_t)last_dim; 
    
    fftw_plan_with_nthreads(n_threads);

    /*
     * Declare variables whose memory is to be freed during cleanup.
     */
    double *k = NULL;
    fftw_plan k_plan  = NULL;
    fftw_plan x1_plan = NULL;
    fftw_plan x2_plan = NULL;
    fftw_plan x3_plan = NULL;
    fftw_plan y1_plan = NULL;
    fftw_plan y2_plan = NULL;
    fftw_plan y3_plan = NULL;
                     
    /*
     * Prepare the kernel
     */
    printf("prepare kernel\n");
    k= malloc(nembed * sizeof(double));
    if (!k) {
        printf("Error: Could not allocate memory for 'k'.\n");
        err = -1;
        goto cleanup;
    }
    double z1, z2, z3, norm;
    ptrdiff_t i1, i2, i3, j1, j2, j3;
#pragma omp parallel for private(z1, z2, z3, norm, i1, i2, i3, j1, j2, j3)
    for (i1 = 0; i1 < dim; i1++) {
        for (i2 = 0; i2 < dim; i2++) {
            for (i3 = 0; i3 < dim; i3++) {
                z1 = (i1 - dim / 2.0 + 0.5) * delta_x;
                z2 = (i2 - dim / 2.0 + 0.5) * delta_x;
                z3 = (i3 - dim / 2.0 + 0.5) * delta_x;
                norm = sqrt(z1 * z1 + z2 * z2 + z3 * z3);
                j1 = WRAP(i1 + 1, dim);
                j2 = WRAP(i2 + 1, dim);
                j3 = WRAP(i3 + 1, dim);
                k[j3 + last_dim * (j2 + dim * j1)] 
                    = z1 / norm * linear_interpolation(norm, k_abscissa[0], k_abscissa[n_k - 1],
                                                       n_k, k_ordinate, 0.0);;
            }
        }
    }

    //print_array_double_3((int[]){dim, dim, last_dim}, k);
    
    /*
     * Plan DFTs and IDFTs
     */
    printf("plan DFTs\n");
    fftw_complex *X1 = (fftw_complex *)x1;
    fftw_complex *X2 = (fftw_complex *)x2;
    fftw_complex *X3 = (fftw_complex *)x3;
    fftw_complex *K  = (fftw_complex *)k;

    k_plan  = fftw_plan_dft_r2c_3d(dim, dim, dim, k,  K,  FFTW_ESTIMATE);
    x1_plan = fftw_plan_dft_r2c_3d(dim, dim, dim, x1, X1, FFTW_ESTIMATE);
    x2_plan = fftw_plan_dft_r2c_3d(dim, dim, dim, x2, X2, FFTW_ESTIMATE);
    x3_plan = fftw_plan_dft_r2c_3d(dim, dim, dim, x3, X3, FFTW_ESTIMATE);
    y1_plan = fftw_plan_dft_c2r_3d(dim, dim, dim, X1, x1, FFTW_ESTIMATE);
    y2_plan = fftw_plan_dft_c2r_3d(dim, dim, dim, X2, x2, FFTW_ESTIMATE);
    y3_plan = fftw_plan_dft_c2r_3d(dim, dim, dim, X3, x3, FFTW_ESTIMATE);

    if (!k_plan || !x1_plan || !x2_plan || !x3_plan || !y1_plan || !y2_plan || !y3_plan) {
        err = AMBIT_FFTW_PLAN_ERROR;
        goto cleanup;
    }

    /*
     * Execute DFTs
     */
    printf("execute DFTs\n");
    fftw_execute(k_plan);
    fftw_execute(x1_plan);
    fftw_execute(x2_plan);
    fftw_execute(x3_plan);

    /*
     * Multiply DFTs
     */
    printf("multiply DFTs\n");
    fftw_complex Y1, Y2, Y3, K1, K2, K3;
    double scale = 1 / ((double)dim * (double)dim * (double)dim);
    ptrdiff_t i;
#pragma omp parallel for private(i, Y1, Y2, Y3, K1, K2, K3, i1, i2, i3)
    for (i1 = 0; i1 < dim; i1++) {
        for (i2 = 0; i2 < dim; i2++) {
            for (i3 = 0; i3 < last_dim / 2; i3++) {
                i  =   i3 + last_dim / 2 * (i2 + dim * i1);
                K1 = K[i3 + last_dim / 2 * (i2 + dim * i1)];
                K2 = K[i3 + last_dim / 3 * (i1 + dim * i2)];
                if (i1 < last_dim / 2) 
                    K3 =      K[i1         + last_dim / 2 * (i2 + dim * i3)];
                else 
                    K3 = conj(K[(dim - i1) + last_dim / 2 * (i2 + dim * i3)]);
                Y1 = scale * (K3 * X2[i] - K2 * X3[i]);
                Y2 = scale * (K1 * X3[i] - K3 * X1[i]);
                Y3 = scale * (K2 * X1[i] - K1 * X2[i]);
                X1[i] = Y1;
                X2[i] = Y2;
                X3[i] = Y3;
            }
        }
    }
    
    /*
     * Execute IDFTs
     */
    printf("execute IDFTs\n");
    fftw_execute(y1_plan);
    fftw_execute(y2_plan);
    fftw_execute(y3_plan);

    /*
     * Cleanup and return
     */
cleanup:
    printf("cleaning up\n");
    free(k);
    fftw_destroy_plan(k_plan);
    fftw_destroy_plan(x1_plan);
    fftw_destroy_plan(x2_plan);
    fftw_destroy_plan(x3_plan);
    fftw_destroy_plan(y1_plan);
    fftw_destroy_plan(y2_plan);
    fftw_destroy_plan(y3_plan);
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
/* int ambit_circular_matrix_convolution_inplace(int rank, */
/*                                               int *kdim, double **k, */
/*                                               int *xdim, double **x,  */
/*                                               int nrows, int ncols, int *matindex) { */
/*     int err = AMBIT_NO_ERROR; */

/*     /\* */
/*      * Inplace DFT of x */
/*      *\/ */
/*     fftw_plan x_plan = NULL; */
/*     fftw_complex **X = (fftw_complex **)malloc(ncols * sizeof(fftw_complex *)); */
/*     for (int i = 0; i < ncols; i++) { */
/*         X[i] = (fftw_complex)x[i]; */
/*         x_plan = fftw_plan_dft_r2c(rank, xdim, x, X, FFTW_ESTIMATE); */
/*         if (!x_plan) { */
/*             err = AMBIT_PLAN_ERROR; */
/*             goto cleanup; */
/*         } */
/*         fftw_execute(x_plan); */
/*         fftw_destroy_plan(x_plan); */
/*     } */
    
/*     /\* */
/*      * Allocate temporary storage */
/*      *\/ */
/*     ptrdiff_t xnembed = 1; */
/*     for (int i = 0; i < rank - 1; i++)  */
/*         xnembed *= xdim[i]; */
/*     xnembed *= 2 * (xdim[rank - 1] / 2 + 1); */
            
/*     double       *t  = NULL; */
/*     fftw_complex *T  = NULL; */
/*     t = (double *)fftw_malloc(xnembed * sizeof(double)); */
/*     T = (fftw_complex *)t; */
    
/*     fftw_complex **Y = (fftw_complex **)malloc(nrows * sizeof(fftw_complex *)); */
/*     for(int i = 0; i < nrows; i++) */
/*         Y[i] = (fftw_complex *)fftw_malloc(xnembed / 2 * sizeof(fftw_complex)); */
        

/*     /\* */
/*      * Cleanup and return */
/*      *\/ */
/* cleanup: */
/*     fftw_destroy_plan(x_plan); */
    


/*     return err; */
/* } */
