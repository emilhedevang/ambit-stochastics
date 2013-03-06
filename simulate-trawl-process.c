#include "ambit-stochastics.h"

int main(int argc, char *argv[]) {
    int err = 0;
    char *usage =
        "Usage: %s {\"N\" mean variance | \"GH\" lambda alpha beta mu delta}"
        " slices.h5 n output.h5\n";
    int argcN  =  7;
    int argcGH = 10;
    
    double lambda, alpha, beta, mu, delta;
    int64_t n = 0;
    char *slice_file  = NULL;
    char *output_file = NULL;
    double *x         = NULL;
    
    void *generator_params;
    int (*generator)(int, gsl_rng **, double, void *, int64_t, double *);

    if (argc != argcN && argc != argcGH) {
        printf(usage, argv[0]);
        err = -1;
        goto cleanup;
    }

    if (strcmp("N", argv[1]) == 0 && argc == argcN) {
        int ctr     = 1;
        mu          = atof(argv[++ctr]);
        delta       = atof(argv[++ctr]);
        slice_file  =      argv[++ctr];
        n           = atol(argv[++ctr]);
        output_file =      argv[++ctr];
        generator_params = (normal_generator_params *)malloc(sizeof(normal_generator_params));
        ((normal_generator_params *)generator_params)->mean     = mu;
        ((normal_generator_params *)generator_params)->variance = delta;
        generator = normal_generator;
    }
    else if (strcmp("GH", argv[1]) == 0 && argc == argcGH) {
        int ctr     = 1;
        lambda      = atof(argv[++ctr]);
        alpha       = atof(argv[++ctr]);
        beta        = atof(argv[++ctr]);
        mu          = atof(argv[++ctr]);
        delta       = atof(argv[++ctr]);
        slice_file  =      argv[++ctr];
        n           = atol(argv[++ctr]);
        output_file =      argv[++ctr];
        generator_params = (univariate_generalised_hyperbolic_generator_params *)
            malloc(sizeof(normal_generator_params));
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->lambda = lambda;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->alpha  = alpha;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->beta   = beta;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->mu     = mu;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->delta  = delta;
        generator = univariate_generalised_hyperbolic_generator;
    }
    else {
        printf("Error: The type of Levy basis must be either normal (N)"
               " or generalised hyperbolic (GH)\n");
        printf(usage, argv[0]);
        err = -1;
        goto cleanup;
    }

    x = calloc(n, sizeof(double));
    if (!x) {
        printf("Error: Could not allocate enough memory.\n");
        err = -1;
        goto cleanup;
    }

    // Read slice volumes
    // Perform simulation
    // Write result
    

  cleanup:
    free(x);
    free(generator_params);
    return err;
}
