#include "ambit-stochastics.h"
#include <argp.h>

#define slice_dataset  "/slice-volumes"
#define output_dataset "/simulation"

//#define DEBUG
#ifdef DEBUG
#define DEBUGPRINT(x)          printf("%s\n", (x))
#else
#define DEBUGPRINT(x)
#endif

struct arguments;
int validate_arguments(struct arguments *p);

const char *argp_program_version     = "simulate-trawl-process v. 0.1";
const char *argp_program_bug_address = "<hedevang@math.au.dk>";

/* Program documentation. */
static char doc[] = "A program to simulate time series from trawl processes";

/* A description of the arguments we accept. */
static char args_doc[] = "";

/* The options we understand. */
static struct argp_option options[] = {
    {0, 0, 0, OPTION_DOC, "Precalculated ambit set:"},
    {"ambit-set",            'A', "file.h5",0, "Use the ambit set specified in 'file.h5'", 0},
    {0, 0, 0, OPTION_DOC, "Parametrised cascade ambit set:"},   
    {"cascade",              'C', 0,        0, "Use the \"cascade\" ambit set",                 0 },
    {"decorrelation-length", 'T', "T",      0, "Decorrelation length of the cascade ambit set", 0 },
    {"cascade-length",       'L', "L",      0, "Cascade length of the cascade ambit set",       0 },
    {"smoothness",           's', "s",      0, "Smoothness parameter of the cascade ambit set", 0 },
    {"resolution",           'r', "dt",     0, "Simulation step size",                          0 },
    {0, 0, 0, OPTION_DOC, "Levy seed:"},
    {"normal",               'N', 0,        0, "Use a normal seed", 0 },
    {"generalised-hyperbolic", 'H', 0,      0, "Use a generalised hyperbolic seed", 0 },
    {"lambda",               'l', "lambda", 0, "Parameter of GH seed", 0 },
    {"alpha",                'a', "alpha",  0, "Parameter of GH seed", 0 },
    {"beta",                 'b', "beta",   0, "Parameter of GH seed", 0 },
    {"mu",                   'm', "mu",     0, "Parameter of GH seed or mean of normal seed", 0 },
    {"delta",                'd', "delta",  0, "Parameter of GH seed or variance of normal seed", 0 },
    {0, 0, 0, OPTION_DOC, "Output:"},
    {"output",               'o', "file.h5",0, "The result is stored in 'file.h5'", 0},
    {"samples",              'n', "n",      0, "Number of samples to generate", 0},
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
    int64_t n;
    char *output_file;
    
    int cascade;
    char *ambit_set_file;
    double T, L, theta, dt;
    
    int normal_seed;
    int generalised_hyperbolic_seed;
    double lambda, alpha, beta, mu, delta;
};


/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
     
    switch (key) {
    case 'A': // ambit-set
        arguments->ambit_set_file = arg;
        break;
    case 'C': // cascade
        arguments->cascade = 1;
        break;
    case 'T': // decorrelation-length
        arguments->T = atof(arg);
        break;
    case 'L': // cascade-length
        arguments->L = atof(arg);
        break;
    case 's': // theta
        arguments->theta = atof(arg);
        break;
    case 'N': // normal seed
        arguments->normal_seed = 1;
        break;
    case 'H': // GH seed
        arguments->generalised_hyperbolic_seed = 1;
        break;
    case 'l': // lambda
        arguments->lambda = atof(arg);
        break;
    case 'a': // alpha
        arguments->alpha = atof(arg);
        break;
    case 'b': // beta
        arguments->beta = atof(arg);
        break;
    case 'm': // mu
        arguments->mu = atof(arg);
        break;
    case 'd': // delta
        arguments->delta = atof(arg);
        break;
    case 'r': // dt
        arguments->dt = atof(arg);
        break;
    case 'n':
        arguments->n = atol(arg);
        break;
    case 'o':
        arguments->output_file = arg;
        break;
    case ARGP_KEY_ARG:
        if (state->arg_num > 0) {
            /* Too many arguments. */
            argp_usage (state);
        }            
        break;
    case ARGP_KEY_END:
        if (validate_arguments(arguments))
            argp_usage (state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/* Initialise argument struct */
void initialise_arguments(struct arguments *p) {
    p->n           = 0;
    p->output_file = NULL;
    
    p->ambit_set_file = NULL;

    p->cascade = 0;
    p->T       = nan("");
    p->L       = nan("");
    p->theta   = nan("");
    p->dt      = nan("");
    
    p->normal_seed                 = 0;
    p->generalised_hyperbolic_seed = 0;
    p->lambda = nan("");
    p->alpha  = nan("");
    p->beta   = nan("");
    p->mu     = nan("");
    p->delta  = nan("");
}

/* Print arguments struct (for debugging) */
void printf_arguments(struct arguments *p) {
    printf("n = %ld\n", p->n);
    printf("output_file = '%s'\n", p->output_file);
    printf("cascade = %d\n", p->cascade);
    printf("ambit_set_file = '%s'\n", p->ambit_set_file);
    printf("(T, L, theta, dt) = (%g, %g, %g, %g)\n", p->T, p->L, p->theta, p->dt);
    printf("normal_seed = %d\n", p->normal_seed);
    printf("generalised_hyperbolic_seed = %d\n", p->generalised_hyperbolic_seed);
    printf("(lambda, alpha, beta, mu, delta) = (%g, %g, %g, %g, %g)\n", p->lambda, p->alpha, p->beta, p->mu, p->delta);
}

/* Validate an argument struct */
int validate_arguments(struct arguments *p) {
    int cascade_OK = 
        p->cascade && 
        isfinite(p->T) && p->T > 0 &&
        isfinite(p->L) && p->L > 1 && 
        isfinite(p->theta) && p->theta > 0 &&
        isfinite(p->dt) && p->dt > 0;
    int ambit_set_file_OK = 
        p->ambit_set_file != NULL;
    int N_OK =
        p->normal_seed &&
        isfinite(p->mu) && 
        isfinite(p->delta) && p->delta > 0;
    int GH_OK = 
        p->generalised_hyperbolic_seed && 
        isfinite(p->lambda) && isfinite(p->alpha) && isfinite(p->beta) && isfinite(p->mu) && isfinite(p->delta) &&
        fabs(p->beta) <= p->alpha && p->delta >= 0;
    
    if ((cascade_OK && ambit_set_file_OK) || (!cascade_OK && !ambit_set_file_OK)) {
        printf("Please specify either a cascade or an ambit set file, but not both.\n");
        return -1;
    }
    if ((N_OK && GH_OK) || (!N_OK && !GH_OK)) {
        printf("Please specify either a normal or a generalised hyperbolic Levy seed, but not both.\n");
        return -1;
    }
    if (p->n < 1) {
        printf("The simulation length must be a positive integer.\n");
        return -1;
    }
    if (!p->output_file) {
        printf("Please specify an output file.\n");
        return -1;
    }
    return 0;
}

double cascade_boundary_function(double t, void *params) {
    double T          = ((double *)params)[0];
    double L          = ((double *)params)[1];
    double smoothness = ((double *)params)[2];
    if (t == 0.0)
        return 1.0;
    else if (0.0 < t && t < T)
        return pow((1 - pow(t / T, smoothness)) / (1 + pow(L * t / T, smoothness)), 1 / smoothness);
    else 
        return 0.0;
}

int calculate_slice_volumes(struct arguments *p, int64_t *n_slices, double **slice_volume) {
    int err = 0;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double T = p->T;
    double L = p->L;
    double theta = p->theta;
    double dt = p->dt;
    double params[3] = {T, L, theta};
    
    *n_slices     = (int64_t)ceil(T / dt);
    *slice_volume = calloc(*n_slices, sizeof(double));
    if (!slice_volume) {
        printf("Error: Could not allocate memory for slice volumes.\n");
        err = -1;
        goto cleanup;
    }
    
    double result   = 0.0;
    double abserror = 0.0;
    gsl_function F;
    F.function = &cascade_boundary_function;
    F.params   = params;
    
    for (int64_t i = 0; i < *n_slices; i++) {
        gsl_integration_qag(&F, i * dt, (i + 1) * dt, 0.0, 1e-7, 1000, GSL_INTEG_GAUSS21, 
                            w, &result, &abserror);
        (*slice_volume)[i] = 2.0 * result;
    }
    
  cleanup:
    gsl_integration_workspace_free(w);
    return err;
}

int read_slice_volumes(struct arguments *p, int64_t *n_slices, double **slice_volume) {
    int     err  = 0;
    hsize_t dim0 = 0;    
    int     slice_rank         = 0;
    hid_t   slice_file_id      = 0;
    hid_t   slice_dataset_id   = 0;
    hid_t   slice_dataspace_id = 0;
    herr_t  status             = 0;

    slice_file_id = H5Fopen(p->ambit_set_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (slice_file_id < 0) {
        printf("Error: Could not open slice file '%s'.\n", p->ambit_set_file);
        err = -1;
        goto cleanup;
    }

    slice_dataset_id = H5Dopen(slice_file_id, slice_dataset, H5P_DEFAULT);
    if (slice_file_id < 0) {
        printf("Error: Could not open slice dataset '%s'.\n", slice_dataset);
        err = -1;
        goto cleanup;
    }

    slice_dataspace_id = H5Dget_space(slice_dataset_id);
    if (slice_dataspace_id < 0) {
        printf("Error: Could not get the dataspace for the dataset '%s'.\n", slice_dataset);
        err = -1;
        goto cleanup;
    }

    slice_rank = H5Sget_simple_extent_ndims(slice_dataspace_id);
    if (slice_rank != 1) {
        printf("Error: The rank of the slice dataset should be one, is %d.\n", slice_rank);
        err = -1;
        goto cleanup;
    }

    H5Sget_simple_extent_dims(slice_dataspace_id, &dim0, NULL);
    *n_slices = (int64_t)dim0;
    *slice_volume = calloc(*n_slices, sizeof(double));
    if (!slice_volume) {
        printf("Error: Could not allocate memory for slices.\n");
        err = -1;
        goto cleanup;
    }
    
    status = H5LTread_dataset_double(slice_file_id, slice_dataset, *slice_volume);
    if (status < 0) {
        printf("Error: Could not read slice dataset '%s'.\n", slice_dataset);
        err = -1;
        goto cleanup;
    }
    
  cleanup:
    if (slice_dataspace_id > 0) H5Sclose(slice_dataspace_id);
    if (slice_dataset_id > 0)   H5Dclose(slice_dataset_id);
    if (slice_file_id > 0)      H5Fclose(slice_file_id);
    return err;
}


int main(int argc, char *argv[]) {
    int err = 0;

    double  *x             = NULL;
    double  *slice_volume  = NULL;
    int64_t  n_slices      = 0;
    herr_t   status        = 0;
    int      n_threads     = 0;
    
    gsl_rng_env_setup();
    const gsl_rng_type  *T   = gsl_rng_default;
    gsl_rng            **rng = NULL;

    void *generator_params = NULL;
    int (*generator)(int, gsl_rng **, double, void *, int64_t, double *);

    hid_t output_file_id     = 0;

    struct arguments *p;
    
    /* Parse command line arguments */
    DEBUGPRINT("parse");
    p = malloc(sizeof(struct arguments));
    if (!p) { printf("Allocation error.\n"); err = -1; goto cleanup; }
    initialise_arguments(p);
    argp_parse (&argp, argc, argv, 0, 0, p);
    
    #ifdef DEBUG
    printf_arguments(p);
    #endif
    
    /* Make seed generators */
    DEBUGPRINT("seed");
    if (p->normal_seed) {
        generator_params = (normal_generator_params *)malloc(sizeof(normal_generator_params));
        ((normal_generator_params *)generator_params)->mean     = p->mu;
        ((normal_generator_params *)generator_params)->variance = p->delta;
        generator = normal_generator;
    }
    else if (p->generalised_hyperbolic_seed) {
        generator_params = (univariate_generalised_hyperbolic_generator_params *)
            malloc(sizeof(univariate_generalised_hyperbolic_generator_params));
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->lambda = p->lambda;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->alpha  = p->alpha;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->beta   = p->beta;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->mu     = p->mu;
        ((univariate_generalised_hyperbolic_generator_params*)generator_params)->delta  = p->delta;
        generator = univariate_generalised_hyperbolic_generator;
    }
    else {
        printf("Error: This should not happen, sorry.\n");
        err = -1;
        goto cleanup;
    }

    /* Get slice volumes */
    DEBUGPRINT("slice");
    if (p->cascade) {
        err = calculate_slice_volumes(p, &n_slices, &slice_volume);
        if (err) {
            printf("Error in calculate_cascade_slice_volumes.\n");
            goto cleanup;
        }
    }
    else if (p->ambit_set_file) {
        err = read_slice_volumes(p, &n_slices, &slice_volume);
        if (err) {
            printf("Error in read_slice_volumes.\n");
            goto cleanup;
        }
    }
    else { 
        printf("Error: This should not happen, sorry.\n"); 
        err = -1; 
        goto cleanup; 
    }

    /* Initialise random number generators */
    DEBUGPRINT("rng");
    n_threads = omp_get_max_threads();
//    printf(" %d\n", n_threads);
    rng       = (gsl_rng**) malloc(n_threads * sizeof(gsl_rng*));
    if (!rng) {
        printf("Error: Could not allocate random number generators.\n");
        err = -1;
        goto cleanup;
    }
    for (int i = 0; i < n_threads; i++) {
//        printf(" %i\n", i);
        rng[i] = gsl_rng_alloc(T);
        gsl_rng_set(rng[i], i + 1);
    }
//    DEBUGPRINT(" ok");
    
    /* Storage for simulation result */
    DEBUGPRINT("x");
    x = calloc(p->n, sizeof(double));
    if (!x) {
        printf("Error: Could not allocate enough memory.\n");
        err = -1;
        goto cleanup;
    }


    /* Perform simulation */
    DEBUGPRINT("sim");
    err = trawl_process(n_threads, rng, generator, generator_params, 
                        n_slices, slice_volume, p->n, x);
    if (err) {
        printf("Error in trawl_process.\n");
        goto cleanup;
    }

    /* Write result */
    DEBUGPRINT("out");
    output_file_id = H5Fcreate(p->output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (output_file_id < 0) {
        printf("Error: could not write to the output file '%s'.\n", p->output_file);
        err = -1;
        goto cleanup;
    }
    status = H5LTmake_dataset_double(output_file_id, output_dataset, 1, (const hsize_t*)&(p->n), x);
    if (status < 0) {
        printf("Error: Could not create the output dataset '%s'.\n", output_dataset);
        err = -1;
        goto cleanup;
    }
    
  cleanup:
    DEBUGPRINT("cleanup");
    if (output_file_id > 0) H5Fclose(output_file_id);
    free(x);
    free(slice_volume);
    free(generator_params);
    for (int i = 0; i < n_threads; i++) 
        gsl_rng_free(rng[i]);
    free(rng);
    free(p);
    return err;
}
