#include "ambit-stochastics.h"

/** Generate multivariate normal variates.
 *
 * n_threads: number of threads to use in OpenMP
 * rng: random number generators, array of length n_threads
 * dim: number of dimensions
 * mean: mean value, array of length dim; zero if NULL
 * covariance: covariance matrix, vector of length dim*dim; identity if NULL
 * n: number of variates to generate
 * x: storage of the result, array of length n*dim
 */

void 
gen_single_mv_normal(const void *param, const void *rng, R *result) {
    R *pp = (R *)param;
    int dim = (int)(pp[0]);
    R *mean = pp + 1;
    R *A    = pp + 1 + dim;
    R *work = pp + 1 + dim + dim * dim;
    
    for (int i = 0; i < dim; ++i) {
        result[i] = mean[i];
        work[i]   = gsl_ran_ugaussian((gsl_rng *)rng);
    }
    cblas_dgemv(CblasColMajor, CblasNoTrans, dim, dim, 1.0, A, dim, work, 1, 1.0, result, 1);
}

int 
generate_multivariate_normal(gsl_rng **rng,
                             struct ambit_dense_array *a,
                             struct ambit_dense_array *mean,
                             struct ambit_dense_array *cov) {
    int err = 0;

    double  *U         = NULL;
    double  *lambda    = NULL;

    if (!ambit_dense_array_valid(a)    || 
        !ambit_dense_array_valid(mean) || 
        !ambit_dense_array_valid(cov)  ||
        !ambit_dense_array_embedded_tightly(mean) ||
        !ambit_dense_array_embedded_tightly(cov)  ||
        mean->rank != 1 || cov->rank != 2 || a->rank > 1) {
        err = -1;
        goto cleanup;
    }
    
    size_t dim = mean->dim[0];
    if (cov->dim[0] != dim || cov->dim[1] != dim || ) {
        err = -1;
        goto cleanup;
    }

    U      = malloc(dim * dim, sizeof(double));
    lambda = malloc(dim      , sizeof(double));
    assert(U && lambda);
        
    err = spectral_decomposition(dim, cov->data, lambda, U);
    if (err) goto cleanup;
    
    for (int i = 0; i < dim; i++)
        if (lambda[i] < 0.0) {
            err = -1;
            goto cleanup;
        }
    
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            U[dim * i + j] *= sqrt(lambda[j]);

    struct gen_single_mv_normal_param *param = malloc(sizeof(struct gen_single_mv_normal_param));
    assert(param);
    
    param->dim = dim;
    param->mean = mean;
    param->A = U;
    
            
#pragma omp parallel private(thread_id, rng0) num_threads(n_threads)
            {
                thread_id = omp_get_thread_num();
                rng0      = rng[thread_id];
                
#pragma omp for
                for (ptrdiff_t j = 0; j < n; j++)
                    x[i + dim * j] = sq * gsl_ran_ugaussian(rng0);
            }
        }
        
        /\* Should be expressed in terms of dgemm when it supports 64 bit matrix sizes *\/
#pragma omp parallel num_threads(n_threads) private(tmp0)
        {
            tmp0 = calloc(dim, sizeof(double));
#pragma omp for
            for (ptrdiff_t j = 0; j < n; j++) {
                cblas_dcopy(dim, x + j * dim, 1, tmp0, 1);
                cblas_dgemv(CblasColMajor, CblasNoTrans, 
                            dim, dim, 1.0, U, dim, tmp0, 1, 0.0, x + j * dim, 1);
            }
            free(tmp0);
        }
    }
    
    
  cleanup:
    free(U);
    free(lambda);
    return err;
}



R
gen_single_uv_normal(const void *param, const void *rng) {
    R *pp    = (R *)param;
    R mean   = pp[0];
    R stddev = pp[1];
    return mean + stddev * gsl_ran_ugaussian((gsl_rng *)rng);
}

int 
generate_univariate_normal(gsl_rng **rng, struct ambit_dense_array *a, R mean, R variance) {
    if (!rng || !ambit_dense_array_valid(a) || variance < 0)
        return -1;
    
    R stddev  = sqrt(variance);
    R param[] = {mean, stddev};
    ambit_dense_array_parallel_fill(a, gen_single_uv_normal,
                                    (const void *)param, sizeof(param),
                                    (const void **)rng, sizeof(gsl_rng *));
    return 0;
}
