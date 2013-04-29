#include "ambit-stochastics.h"

int
ambit_dft_r2c_inplace(struct ambit_dense_array *a, int sign) {
    int       err  = 0;
    fftw_plan plan = NULL;
    
    if (!ambit_dense_array_fftw_r2c_embedded(a)) {
        fprintf(stderr, "%s: Invalid array or inappropriate embedding.\n", __func__);
        err = -1;
        goto cleanup;
    }
    
    switch (sign) {
    case FFTW_FORWARD:
        plan = fftw_plan_many_dft_r2c(a->rank, (const int *)a->dim, 1, 
                                      a->data, (const int *)a->dimembed, 1, 0, 
                                      (C *)a->data, (const int *)a->dimembed, 1, 0, 
                                      FFTW_ESTIMATE);
        break;
    case FFTW_BACKWARD:
        plan = fftw_plan_many_dft_c2r(a->rank, (const int *)a->dim, 1, 
                                      (C *)a->data, (const int *)a->dimembed, 1, 0, 
                                      a->data, (const int *)a->dimembed, 1, 0, 
                                      FFTW_ESTIMATE);
        break;
    default:
        fprintf(stderr, "%s: Sign must be FFTW_FORWARD or FFTW_BACKWARD.\n", __func__);
        err = -1;
        goto cleanup;
    }

    if (!plan) {
        fprintf(stderr, "%s: DFT planning failed.\n", __func__);
        err = -1;
        goto cleanup;
    }
    
    fftw_execute(plan);

  cleanup:
    fftw_destroy_plan(plan);
    return err;
}

