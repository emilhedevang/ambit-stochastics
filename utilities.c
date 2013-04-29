#include "ambit-stochastics.h"


gsl_rng **
ambit_rng_alloc(int n) {
    assert(n>0);
    
    const gsl_rng_type *rng_type = NULL;
    gsl_rng **rng = NULL;
    
    gsl_rng_env_setup();
    
    rng_type = gsl_rng_default;
    rng = malloc((n + 1) * sizeof(gsl_rng *));
    assert(rng);
    
    for (int i = 0; i < n; ++i) {
        rng[i] = gsl_rng_alloc(rng_type);
        assert(rng[i]);
        gsl_rng_set(rng[i], gsl_rng_default_seed + i);
    }
    rng[n] = NULL;
    return rng;
}

void 
ambit_rng_free(gsl_rng **rng) {
    while(*rng)
        gsl_rng_free(*rng++);
}


void print_array_double_3(int *dims, double *x) {
    for (int i0 = 0; i0 < dims[0]; i0++) {
        for (int i1 = 0; i1 < dims[1]; i1++) {
            for (int i2 = 0; i2 < dims[2]; i2++)
                printf(" %8.5f", x[i2 + dims[2] * (i1 + dims[1] * i0)]);
            printf("\n");
        }
        printf("\n");
    }
}
    

double linear_interpolation(double x, double a, double b, int n, double *tbl, double y0) {
    double delta = (b - a) / (n - 1);
    if (x < a || x > b)
        return y0;
    if (x == b)
        return tbl[n - 1];
    int i = (int)floor((x - a) / delta);
    double t = x - a - i * delta;
    return (1 - t) * tbl[i] + t * tbl[i + 1];
}

int spectral_decomposition(int dim, double *a, double *lambda, double *u) {
    int err              = 0;
    lapack_int  info     = 0;
    lapack_int  ev_found = 0;
    double     *a_copy   = NULL;
    lapack_int *isuppz   = NULL;
    
    if (dim < 1 || !a || !lambda || !u) {
        err = -1;
        goto cleanup;
    }
    
    a_copy = calloc(dim * dim, sizeof(double));
    isuppz = calloc(2 * dim  , sizeof(lapack_int));
    
    if (!a_copy || !isuppz) {
        err = -1;
        goto cleanup;
    }
    
    cblas_dcopy(dim * dim, a, 1, a_copy, 1);
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'L', 
                          dim, a_copy, dim, 
                          0.0, 0.0, 0, 0, 0.0, 
                          &ev_found, lambda, u, dim, isuppz);
    
    if (info || ev_found < dim) {
        err = -1;
        goto cleanup;
    }
    
  cleanup:
    free(a_copy);
    free(isuppz);
    return err;
}
