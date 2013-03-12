#include "ambit-stochastics.h"

int main(int argc, char *argv[]) {
    
    int err = 0;
    
    int64_t  n         = 0;
    double  *x         = NULL;
    int      n_threads = 0;

    int output = false;

    gsl_rng_env_setup();
    const gsl_rng_type  *T   = gsl_rng_default;
    gsl_rng            **rng = NULL;

    if (argc != 6) {
        err = -1;
        printf("Usage: %s n s1 s2 s3 s4\n", argv[0]);
        goto cleanup;
    }
    
    double slice_volumes[4];
    int ctr = 0;
    n         = atol(argv[++ctr]);
    output    = n > 0;
    n         = labs(n);
    for (int i = 0; i < 4; i++)
        slice_volumes[i] = atof(argv[++ctr]);
    
    x         = calloc(n, sizeof(double));
    if (!x) {err = -1; goto cleanup;}
    
    n_threads = omp_get_max_threads();
    rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
    for (int i = 0; i < n_threads; i++) {
        rng[i] = gsl_rng_alloc(T);
        gsl_rng_set(rng[i], i + 1);
    }
    
    err = trawl_process(0, NULL, constant_generator, NULL, 4, slice_volumes, n, x);
    if (err) {
        printf("Error in trawl_process.\n");
        goto cleanup;
    }
    
    if (output)
        for (int64_t j = 0; j < n; j++)
            printf(" %f\n", x[j]);
    
  cleanup:
    free(x);
    return err;

    
}
