#include "ambit-stochastics.h"
#include <argp.h>


#define DEBUG
#ifdef DEBUG
#define DEBUGPRINT(x)  printf("%s\n", (x))
#else
#define DEBUGPRINT(x)
#endif


struct arguments;
int valid_arguments(struct arguments *p);

/* 
 * Argp setup
 */

const char *argp_program_version     = "simulate-homogeneous-levy-basis v. 0.1";
const char *argp_program_bug_address = "<hedevang@math.au.dk>";

/* Program documentation. */
static char doc[] = "A program to simulate homogeneous Levy bases";

/* A description of the arguments we accept. */
static char args_doc[] = "input.h5 output.h5";

/* The options we understand. */
static struct argp_option options[] = {
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
    char *input_file;
    char *output_file;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *p = state->input;
    switch (key) {
    case ARGP_KEY_ARG:
        switch (state->arg_num) {
        case 0: p->input_file = arg;  break;
        case 1: p->output_file = arg; break;
        default: argp_usage(state);
        }
        break;
    case ARGP_KEY_END:
        if (!valid_arguments(p))
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
    p->input_file  = NULL;
    p->output_file = NULL;
}

/* Print arguments struct (for debugging) */
void print_arguments(struct arguments *p) {
    printf("input_file = \"%s\"\n", p->input_file);
    printf("output_file = \"%s\"\n", p->output_file);
}

/* Validate an argument struct */
int valid_arguments(struct arguments *p) {
    if (p->input_file && p->output_file)
        return true;
    else 
        return false;
}

struct input_config {
    char   *rng_family;
    int     rng_seed;
    char   *levy_seed_family;
    int     levy_seed_dimension;
    void   *levy_seed_parameters;
    int     levy_basis_rank;
    int    *levy_basis_dimension;
    double *levy_basis_resolution;
};


struct normal_seed_parameters {
    int     dimension;
    double *mean;
    double *covariance;
};

void print_input_config(struct input_config *p) {
    printf("rng_family = \"%s\"\n", p->rng_family);
    printf("rng_seed  = %i\n", p->rng_seed);
    printf("levy_seed_family = \"%s\"\n", p->levy_seed_family);
    printf("levy_seed_dimension = %i\n", p->levy_seed_dimension);
    if (strcmp(p->levy_seed_family, "normal") == 0) {
        struct normal_seed_parameters *q = p->levy_seed_parameters;
        printf("normal_seed_mean =");
        for (int i = 0; i < q->dimension; i++) 
            printf(" %g", q->mean[i]);
        printf("\n");
        printf("normal_seed_covariance =");
        for (int i = 0; i < q->dimension * q->dimension; i++) 
            printf(" %g", q->covariance[i]);
        printf("\n");
    }
    printf("levy_basis_rank = %i\n", p->levy_basis_rank);
    printf("levy_basis_dimension =");
    for (int i = 0; i < p->levy_basis_rank; i++)
        printf(" %i", p->levy_basis_dimension[i]);
    printf("\n");
    printf("levy_basis_resolution =");
    for (int i = 0; i < p->levy_basis_rank; i++)
        printf(" %g", p->levy_basis_resolution[i]);
    printf("\n");
}


herr_t read_normal_seed_parameters(hid_t loc_id, struct normal_seed_parameters *p) {
    herr_t err = 0;
    err = H5LTread_dataset_int(loc_id, "/levy_seed_dimension", &p->dimension);
    p->mean       = malloc(p->dimension * sizeof(double));
    p->covariance = malloc(p->dimension * p->dimension * sizeof(double));
    err = H5LTread_dataset_double(loc_id, "/levy_seed_parameters/mean",       p->mean);
    err = H5LTread_dataset_double(loc_id, "/levy_seed_parameters/covariance", p->covariance);
    return err;
}

int read_input_config(char *input_file, struct input_config *cfg) {
    herr_t   err = 0;
    size_t   type_size;
    hid_t    input_file_id;
    void    *p = NULL;

    cfg->rng_family            = NULL;
    cfg->levy_seed_family      = NULL;
    cfg->levy_seed_parameters  = NULL;
    cfg->levy_basis_dimension  = NULL;
    cfg->levy_basis_resolution = NULL;
    
    input_file_id = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);
        
    err = H5LTget_dataset_info(input_file_id, "/rng_family", NULL, NULL, &type_size);
    cfg->rng_family = malloc(type_size * sizeof(char));
    err = H5LTread_dataset_string(input_file_id, "/rng_family", cfg->rng_family);
    
    err = H5LTread_dataset_int(input_file_id, "/rng_seed", &cfg->rng_seed);
    
    err = H5LTget_dataset_info(input_file_id, "/levy_seed_family", NULL, NULL, &type_size);
    cfg->levy_seed_family = malloc(type_size * sizeof(char));
    err = H5LTread_dataset_string(input_file_id, "/levy_seed_family", cfg->levy_seed_family);

    err = H5LTread_dataset_int(input_file_id, "/levy_seed_dimension", &cfg->levy_seed_dimension);
    
    if (strcmp(cfg->levy_seed_family, "normal") == 0) {
        p = (struct normal_seed_parameters *)malloc(sizeof(struct normal_seed_parameters));
        err = read_normal_seed_parameters(input_file_id, p);
        cfg->levy_seed_parameters = (void *)p;
    }
    else if (strcmp(cfg->levy_seed_family, "generalised hyperbolic") == 0) {
        printf("Error: GH not yet supported.\n");
        err = -1;
        goto cleanup;
    }
    else {
        printf("Error: Unsupported Levy seed.\n");
        err = -1;
        goto cleanup;
    }
    
    err = H5LTread_dataset_int(input_file_id, "/levy_basis_rank", &cfg->levy_basis_rank);
    cfg->levy_basis_dimension = malloc(cfg->levy_basis_rank * sizeof(int));
    cfg->levy_basis_resolution = malloc(cfg->levy_basis_rank * sizeof(double));
    err = H5LTread_dataset_int(input_file_id, "/levy_basis_dimension", cfg->levy_basis_dimension);
    err = H5LTread_dataset_double(input_file_id, "/levy_basis_resolution", cfg->levy_basis_resolution);
    
  cleanup:
    if (input_file_id >= 0) H5Fclose(input_file_id);
    return (int)err;
}

int write_output(char *output_file, struct input_config *p, double *x) {
    herr_t err = 0;
    hid_t output_file_id;
    int rank;
    hsize_t *dims = NULL;

    output_file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (output_file_id <= 0) {
        printf("Error: Could not create \"%s\".\n", output_file);
        err = -1;
        goto cleanup;
    }
    rank = p->levy_basis_rank + 1;
    dims = malloc(rank * sizeof(hsize_t));
    for (int i = 0; i < rank - 1; i++)
        dims[i] = p->levy_basis_dimension[i];
    dims[rank - 1] = p->levy_seed_dimension;
    err = H5LTmake_dataset_double(output_file_id, "/levy_basis_realization", rank, dims, x);

  cleanup:
    return (int)err;
}


int main(int argc, char *argv[]) {
    int err = 0;
    double *m = NULL;
    double *c = NULL;
    double *x = NULL;

    gsl_rng_env_setup();
    const gsl_rng_type  *T   = gsl_rng_default;
    gsl_rng            **rng = NULL;

    DEBUGPRINT("### Parsing arguments");
    struct arguments args;
    initialise_arguments(&args);
    argp_parse (&argp, argc, argv, 0, 0, &args);
#ifdef DEBUG
    print_arguments(&args);
#endif

    DEBUGPRINT("### Reading input configuration");
    struct input_config cfg;
    read_input_config(args.input_file, &cfg);
#ifdef DEBUG
    print_input_config(&cfg);
#endif
    
    DEBUGPRINT("### Generating realization");
    int  n_threads = omp_get_max_threads();
    
    rng = malloc(n_threads * sizeof(gsl_rng *));
    if (!rng) {
        printf("Error: Could not allocate memory for 'rng'.\n");
        err = -1;
        goto cleanup;
    }
    for (int i = 0; i < n_threads; i++) {
        rng[i] = gsl_rng_alloc(T);
        gsl_rng_set(rng[i], cfg.rng_seed + i);
    }

    ptrdiff_t  n   = 1;
    double     vol = 1.0;
    for (int i = 0; i < cfg.levy_basis_rank; i++) {
        n *= (ptrdiff_t)cfg.levy_basis_dimension[i];
        vol *= cfg.levy_basis_resolution[i];
    }
    
    x = malloc(n * cfg.levy_seed_dimension * sizeof(double));
    if (!x) {
        printf("Error: Could not allocate memory for 'x'.\n");
        err = -1;
        goto cleanup;
    }
    
    if (strcmp(cfg.levy_seed_family, "normal") == 0) {
        struct normal_seed_parameters *q = cfg.levy_seed_parameters;
        int dim = cfg.levy_seed_dimension;
        m = calloc(dim, sizeof(double));
        c = calloc(dim * dim, sizeof(double));
        for (int i = 0; i < dim; i++)
            m[i] = q->mean[i] * vol;
        for (int i = 0; i < dim * dim; i++)
            c[i] = q->covariance[i] * vol;
        err = generate_multivariate_normal(n_threads, rng, dim, m, c, n, x);
        if (err) {
            printf("Error in generate_multivariate_normal.\n");
            goto cleanup;
        }
        
    }
    else {
        printf("Error: The Levy seed family \"%s\" is not yet supported.\n", cfg.levy_seed_family);
        err = -1;
        goto cleanup;
    }

    DEBUGPRINT("### Writing output");
    write_output(args.output_file, &cfg, x);
    
  cleanup:
    free(c);
    free(m);
    free(x);
    for (int i = 0; i < n_threads; i++)
        gsl_rng_free(rng[i]);
    free(rng);
    return err;
}
