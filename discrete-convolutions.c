#include "ambit-stochastics.h"

/** 
 * Calculates the circular convolution of `k` and `x`, overwriting `x`
 * with the result.
 *
 * The array `x` is assumed to be padded along the last dimension with
 * 2 extra elements if `dim[rank - 1]` is even and 1 extra element if
 * `dim[rank - 1]` is odd. For example, 
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

int 
ambit_circular_convolution_inplace(struct ambit_dense_array *k, struct ambit_dense_array *x) {
    int  err = 0;
    struct ambit_dense_array *y      = NULL;
    size_t                   *offset = NULL;
    
    //
    // Initial checks
    //
    if (!ambit_dense_array_dims_leq(k, x)) {
        fprintf(stderr, 
                "%s: Invalid arrays or dimensions of 'k' exceed those of 'x'.\n", __func__);
        err = -1;
        goto cleanup;
    }

    if(!ambit_dense_array_fftw_r2c_embedded_tightly(x)) {
        fprintf(stderr, 
                "%s: Embedding of 'x' is not tightly the FFTW R2C format.\n", __func__);
        err = -1;
        goto cleanup;
    }

    y      = ambit_dense_array_malloc(x->rank, x->dim, x->dimembed);
    offset = malloc(x->rank * sizeof(size_t));
    if (!ambit_dense_array_valid(y) || !offset) {
        fprintf(stderr,
                "%s: Could not allocate temporary memory.\n", __func__);
        err = -1;
        goto cleanup;
    }

    //
    // Circular embedding of kernel
    //
    for (int i = 0; i < x->rank; ++i)
        offset[i] = x->dim[i] - k->dim[i] + 1;
    ambit_dense_array_circular_embed(k, y, offset);

    //
    // DFTs
    //
    err = ambit_dft_r2c_inplace(y, FFTW_FORWARD);
    if (err) {
        fprintf(stderr, "%s: DFT of 'k' failed.\n", __func__);
        goto cleanup;
    }
    err = ambit_dft_r2c_inplace(x, FFTW_FORWARD);
    if (err) {
        fprintf(stderr, "%s: DFT of 'x' failed.\n", __func__);
        goto cleanup;
    }

    //
    // Multiply DFTs
    //
    C *Y = (C *)y->data;
    C *X = (C *)x->data;
    size_t n     = x->nembed / 2;
    double scale = 1 / (double)x->n;
    for (ptrdiff_t i = 0; i < n; ++i) 
        X[i] *= scale * Y[i];

    //
    // Backward DFT
    //
    err = ambit_dft_r2c_inplace(x, FFTW_BACKWARD);
    if (err) {
        fprintf(stderr, "%s: IDFT of 'x' failed.\n", __func__);
        goto cleanup;
    }

  cleanup:
    ambit_dense_array_free(y);
    free(offset);
    return err;
}
    
