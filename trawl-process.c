#include "ambit-stochastics.h"

#define SIGN(x) ((int)(((x) < 0) ? -1 : ((x) > 0)))


/*
 * Constant generator for debugging 
 */

int constant_generator(int n_threads, gsl_rng **rng, 
                       double volume, void *p, 
                       int64_t n, double *x) {
    for (int64_t j = 0; j < n; j++)
        x[j] = j * volume;
    return 0;
}
    
/* 
 * Normal generator
 */

int normal_generator(int n_threads, gsl_rng **rng, 
                     double volume, void *p, 
                     int64_t n, double *x) {
    if (!p) return -1;
    normal_generator_params *params = (normal_generator_params *)p;
    return generate_univariate_normal(n_threads, rng,
                                      params->mean * volume, params->variance * volume,
                                      n, x);
}

/* 
 * Generalised hyperbolic generator
 */

int univariate_generalised_hyperbolic_generator(int n_threads, gsl_rng **rng, 
                                                double volume, void *p, 
                                                int64_t n, double *x) {
    univariate_generalised_hyperbolic_generator_params *params =
        (univariate_generalised_hyperbolic_generator_params *)p;
    return generate_univariate_generalised_hyperbolic(
        n_threads, rng,
        params->lambda, params->alpha, params->beta, params->mu * volume, params->delta * volume,
        n, x);
}

/*
 * Trawl processes
 */

int trawl_process(
    int n_threads, gsl_rng **rng,
    int (*generator)(int, gsl_rng **, double, void *, int64_t, double *),
    void *generator_params,
    int64_t n_slices, double *slice_volume,
    int64_t n, double *x) {
   
    int     err = 0;
    double *y   = NULL;
    double *v   = NULL;
    
    y = calloc(n + n_slices - 1, sizeof(double));
    v = calloc(n_slices, sizeof(double));
    if (!y || !v) {err = -1; goto cleanup;}
    
    /* Test for monotonicity of slice volumes */
    int64_t n_pos = 0, n_neg = 0;
    for (int64_t i = 1; i < n_slices; i++)
        switch (SIGN(slice_volume[i] - slice_volume[i - 1])) {
        case +1: n_pos++; break;
        case -1: n_neg++; break;
        case  0: break;
        }
    if (n_pos > 0 && n_neg > 0) {err = -1; goto cleanup;}
    int monotonicity = (n_neg > 0 ? -1 : (n_pos > 0));
    
    /* Calulate volumes of pieces */
    switch (monotonicity) {
    case 0:
        v[0] = slice_volume[0];
        break;
    case 1:
        v[0] = slice_volume[0];
        for (int64_t i = 1; i < n_slices; i++)
            v[i] = slice_volume[i] - slice_volume[i - 1];
        break;
    case -1:
        v[0] = slice_volume[n_slices - 1];
        for (int64_t i = 1; i < n_slices; i++)
            v[i] = slice_volume[n_slices - 1 - i] - slice_volume[n_slices - 1 - (i - 1)];
        break;
    }
        
    double  acc    = 0.0;
    int64_t offset = 0;
    
    for (int64_t i = 0; i < n_slices; i++) {
        err = (*generator)(n_threads, rng, v[i], generator_params, n + n_slices - i - 1, y);
        if (err) goto cleanup;
        acc    = 0.0;
        offset = monotonicity == 1 ? i : 0;
        for (int64_t j = 0; j < n_slices - i; j++)
            acc += y[offset + j];
        for (int64_t j = 0; j < n; j++) {
            x[j] += acc;
            acc  += y[offset + n_slices - i + j] - y[offset + j];
        }
    }

  cleanup:
    free(y);
    free(v);
    return err;
}
