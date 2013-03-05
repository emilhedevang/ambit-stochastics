#define _POSIX_C_SOURCE 200809L

#include "discrete-convolutions.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <time.h>

double elapsed(struct timespec tic, struct timespec toc) {
    return 
        (double)(toc.tv_sec - tic.tv_sec) 
        + (double)(toc.tv_nsec - tic.tv_nsec) / (double)1000000000;
}

void test_increment_subscript() {
    int rank = 3;
    int dim[3] = {3,4,5};
    int sub[3] = {0,0,0};
    int n = 1;
    for (int i = 0; i < rank; i++)
        n *= dim[i];

    for (int i = 0; i < n; i++) {
        printf("%d: %d %d %d\n", i, sub[0], sub[1], sub[2]);
        increment_subscript(rank, dim, sub);
    }
}    

void test_embed_array_double() {
    int rank = 2;
    int xdim[2] = {2,3};
    int ydim[2] = {7,8};
    int offset[2] = {3,2};
    int xn = xdim[0] * xdim[1];
    int yn = ydim[0] * ydim[1];
    double *x = (double *)malloc(xn * sizeof(double));
    double *y = (double *)malloc(yn * sizeof(double));
    for (int i = 0; i < xn; i++) 
        x[i] = (double)(i+1);
    for (int i = 0; i < yn; i++)
        y[i] = 0.0;
    embed_array_double(rank, xdim, xdim, x, ydim, ydim, y, offset);
    for (int i0 = 0; i0 < ydim[0]; i0++) {
        for (int i1 = 0; i1 < ydim[1]; i1++) {
            printf(" %g", y[i1 + i0 * ydim[1]]);
        }
        printf("\n");
    }
}

void test_ambit_circular_convolution_inplace() {
    int rank = 2;
    int kdim[2] = {4, 5};
    int x0dim[2] = {10, 12};
    int x1dim[2] = {10, 14};
    double k[]  = {-9, -9, 4, 2, -2, -1, -4, 2, 4, -1, -7, 4, -8, 6, 3, 7, 7, -6, 4, 8};
    double x0[]  = {-4, 5, -9, -4, 9, -3, 3, 6, -5, -3, -2, -9, -5, 0, -1, -3, 3, 9, -2, 
                   -5, -1, -9, -9, 9, 3, -6, -5, -1, 2, 0, -9, -3, -3, 8, 6, -4, 5, -1, 
                   3, 6, -5, -2, -4, 1, 5, -4, 5, 5, -7, -1, 8, -5, 3, -4, -9, 3, 8, 3, 
                   6, 9, -2, 1, -7, -4, -4, -2, 5, 1, 2, 3, 0, -4, -2, 0, -9, 6, 6, -8, 
                   4, 6, -4, -1, -8, -3, -4, 4, 7, -1, -5, 8, 7, 4, 9, 1, -2, 0, 1, 7, 
                   4, -8, -2, 9, 6, 2, 4, -4, 6, 4, -5, 0, 1, 5, -3, 3, 1, -8, -3, 6, 2,
                   9};
    double kx[] = {-7, 108, -83, 94, 46, -37, 16, -190, -73, -140, 33, -255, 1, 120, 
                   163, -82, -52, -199, -183, -156, -48, 78, -91, -3, 167, -50, -48, 
                   -123, -58, -61, -73, 32, -56, 151, -33, 38, -235, 49, 129, -162, 1, 
                   -4, -94, 208, 190, -74, 178, 248, -2, -27, -219, -148, -48, -29, 35, 
                   205, 86, -71, 110, 18, 183, -127, -353, 39, -41, -31, 200, -99, -8, 
                   -27, -137, -199, 174, -86, -196, 144, 89, -112, -74, -76, 19, 9, 
                   -169, -150, -178, 50, 207, -171, 52, 328, -83, 162, 122, -56, 148, 
                   145, -21, -102, 144, 90, 40, 207, 212, -62, 38, 31, 200, -17, -56, 
                   -27, -24, 157, 13, -30, 36, 56, -38, -15, 271, 17};

    double *x1 = (double *)malloc(x1dim[0] * x1dim[1] * sizeof(double));
    int offset[2] = {0, 0};
    embed_array_double(rank, x0dim, x0dim, x0, x1dim, x1dim, x1, offset);
    int err = ambit_circular_convolution_inplace(rank, kdim, k, x0dim, x1);
    printf("err = %d\n", err);
    double max = 0;
    for (int i0 = 0; i0 < x0dim[0]; i0++) {
        for (int i1 = 0; i1 < x0dim[1]; i1++) {
            int sub[2] = {i0, i1};
            max = fmax(max, fabs(kx[linear_index(rank, x0dim, x0dim, sub)]
                                 - x1[linear_index(rank, x0dim, x1dim, sub)]));
        }
    };
    printf("max = %g\n", max);

}

void test_3dim_circular_convolution(int argc, char* argv[]) {
    int kdim0, xdim0;
    if (argc == 2 + 1) {
        kdim0 = atoi(argv[1]);
        xdim0 = atoi(argv[2]);
    }
    else {
        printf("Error: You must supply three arguments.\n");
        return;
    }
    int kdim[3] = {kdim0, kdim0, kdim0};
    int xdim[3] = {xdim0, xdim0, xdim0};

    
    struct timespec tic, toc;
    int err;
    int n_threads = omp_get_max_threads();
    
    printf("Size of k: %d * %d * %d\n", kdim[0], kdim[1], kdim[2]);
    printf("Size of x: %d * %d * %d\n", xdim[0], xdim[1], xdim[2]);
    printf("Number of threads: %i\n", n_threads);

    err = fftw_init_threads();
    if (err == 0) {
        printf("Error in 'fftw_init_threads'.\n");
        return;
    }
    
    /* 
       Setup random number generators
    */
    const gsl_rng_type *T;
    gsl_rng **r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = malloc(n_threads * sizeof(gsl_rng*));
    for (int i = 0; i < n_threads; i++) { 
        r[i] = gsl_rng_alloc (T);
        gsl_rng_set(r[i], i + 1);
    }
    
    /* 
       Arrays
    */
    ptrdiff_t n_k = kdim[0] * kdim[1] * kdim[2];
    ptrdiff_t n_x = xdim[0] * xdim[1] * 2 * (xdim[2] / 2 + 1);
    double *k = NULL;
    double *x = NULL;
    k = fftw_alloc_real(n_k);
    x = fftw_alloc_real(n_x);

    if (!k || !x) {
        printf("Error: Not enough memory.\n");
        return;
    }
    
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_k; i++)
        k[i] = 0.0;
    
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n_x; i++)
        x[i] = 0.0;
    
    clock_gettime(CLOCK_REALTIME, &tic);
    int thread_id;
    gsl_rng *r0;
#pragma omp parallel shared(k, x)  private(thread_id, r0) num_threads(n_threads)
    {
        thread_id = omp_get_thread_num();
        r0 = r[thread_id];
        
#pragma omp for
        for (ptrdiff_t i = 0; i < n_k; i++)
            k[i] = gsl_ran_ugaussian(r0);
        
#pragma omp for
        for (ptrdiff_t i = 0; i < n_x; i++)
            x[i] = gsl_ran_ugaussian(r0);
    }
    clock_gettime(CLOCK_REALTIME, &toc);
    printf("Filling k and x with random numbers took %g seconds\n", elapsed(tic, toc));
    
    clock_gettime(CLOCK_REALTIME, &tic);
    err = ambit_circular_convolution_inplace(3, kdim, k, xdim, x);
    clock_gettime(CLOCK_REALTIME, &toc);
    if (!err)
        printf("Convolving k and x took %g seconds\n", elapsed(tic, toc));
    else
        printf("Some error occurred while convolving k and x: err = %d\n", err);

    fftw_cleanup_threads();
}


int main (int argc, char *argv[]) {
    
//    test_increment_subscript();
//    test_embed_array_double();
//    test_ambit_circular_convolution_inplace();
    test_3dim_circular_convolution(argc, argv);
    
    return 0;

}
