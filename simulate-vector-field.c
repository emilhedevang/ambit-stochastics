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

const char *argp_program_version     = "simulate-vector-field v. 0.1";
const char *argp_program_bug_address = "<hedevang@math.au.dk>";

/* Program documentation. */
static char doc[] = "A program to simulate vector fields";

/* A description of the arguments we accept. */
static char args_doc[] = "kernel.h5 levy-basis.h5 output.h5";

/* The options we understand. */
static struct argp_option options[] = {
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
    char *kernel_file;
    char *levy_basis_file;
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
        case 0: p->kernel_file     = arg; break;
        case 1: p->levy_basis_file = arg; break;
        case 2: p->output_file     = arg; break;
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
    p->kernel_file     = NULL;
    p->levy_basis_file = NULL;
    p->output_file     = NULL;
}

/* Print arguments struct (for debugging) */
void print_arguments(struct arguments *p) {
    printf("kernel_file = \"%s\"\n",     p->kernel_file);
    printf("levy_basis_file = \"%s\"\n", p->levy_basis_file);
    printf("output_file = \"%s\"\n",     p->output_file);
}

/* Validate an argument struct */
int valid_arguments(struct arguments *p) {
    if (p->kernel_file && p->levy_basis_file && p->output_file)
        return true;
    else
        return false;
}

int main(int argc, char *argv[]) {
    herr_t err = 0;
    
    int n_threads = omp_get_max_threads();

    hid_t kernel_file_id = 0;
    hid_t levy_basis_file_id = 0;
    hid_t levy_basis_dataset_id = 0;
    hid_t levy_basis_dataspace_id = 0;
    hid_t output_file_id      = 0;
    hid_t output_dataset_id   = 0;
    hid_t output_dataspace_id = 0;
    hid_t memspace = 0;

    hsize_t n_k = 0;
    double *tmp = NULL;
    double *k_abscissa = NULL;
    double *k_ordinate = NULL;
    double *x1 = NULL;
    double *x2 = NULL;
    double *x3 = NULL;
    
    DEBUGPRINT("### Parsing arguments");
    struct arguments args;
    initialise_arguments(&args);
    argp_parse (&argp, argc, argv, 0, 0, &args);
#ifdef DEBUG
    print_arguments(&args);
#endif
    
    DEBUGPRINT("### Reading kernel");
    kernel_file_id = H5Fopen(args.kernel_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (kernel_file_id <= 0) {
        printf("Error: Could not open \"%s\".\n", args.kernel_file);
        err = -1;
        goto cleanup;
    }
    
    err = H5LTget_dataset_info(kernel_file_id, "/abscissa", &n_k, NULL, NULL);
    if (err < 0) {
        printf("Error: Could not read dataset info.\n");
        goto cleanup;
    }
    printf("n_k = %i\n", (int)n_k);
    k_abscissa = malloc(n_k * sizeof(double));
    k_ordinate = malloc(n_k * sizeof(double));
    err = H5LTread_dataset_double(kernel_file_id, "/abscissa", k_abscissa);
    err = H5LTread_dataset_double(kernel_file_id, "/ordinate", k_ordinate);
    
    DEBUGPRINT("### Reading Levy basis");
    hsize_t dims[4];
    hsize_t offset[4];
    hsize_t count[4];
    levy_basis_file_id      = H5Fopen(args.levy_basis_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    levy_basis_dataset_id   = H5Dopen(levy_basis_file_id, "/levy_basis_realization", H5P_DEFAULT);
    levy_basis_dataspace_id = H5Dget_space(levy_basis_dataset_id);
    err = H5Sget_simple_extent_dims(levy_basis_dataspace_id, dims, NULL);
    
    if (dims[0] != dims[1] || dims[1] != dims[2]) {
        printf("Error: The three dimensions must be equal.\n");
        err = -1;
        goto cleanup;
    }

    hsize_t dims_pad[3];
    dims_pad[0] = dims[0];
    dims_pad[1] = dims[1];
    dims_pad[2] = 2 * (dims[2] / 2 + 1);

    hsize_t n_x = dims_pad[0] * dims_pad[1] * dims_pad[2];
    x1 = malloc(n_x * sizeof(double));
    x2 = malloc(n_x * sizeof(double));
    x3 = malloc(n_x * sizeof(double));

    double *x[] = {x1, x2, x3};
    if (!x1 || !x2 || !x3) {
        printf("Error: Could not allocate memory for the Levy basis.\n");
        err = -1;
        goto cleanup;
    }

#pragma omp parallel for 
    for (ptrdiff_t i = 0; i < n_x; i++) {
        x1[i] = 0.0;
        x2[i] = 0.0;
        x3[i] = 0.0;
    }    

    /* Define memory dataspace */
    memspace = H5Screate_simple(3, dims_pad, NULL);
    offset[0] = offset[1] = offset[2] = 0;
    count[0] = dims[0];
    count[1] = dims[1];
    count[2] = dims[2];
    /* Define hyperslab in the memory dataspace */
    err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    for (int j = 0; j < 3; j++) {
        /* Define hyperslap in the file dataspace */
        offset[3] = j;
        count[3] = 1;
        err = H5Sselect_hyperslab(levy_basis_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        /* Read data from hyperslab */
        err = H5Dread(levy_basis_dataset_id, H5T_NATIVE_DOUBLE, 
                      memspace, levy_basis_dataspace_id,
                      H5P_DEFAULT, x[j]);
        if (err < 0) {
            printf("Error: Could not read hyperslab.\n");
            err = -1;
            goto cleanup;
        }
        
        /* for (int i0 = 0; i0 < dims_pad[0]; i0++) { */
        /*     for (int i1 = 0; i1 < dims_pad[1]; i1++) { */
        /*         for (int i2 = 0; i2 < dims_pad[2]; i2++) */
        /*             printf(" %g", x[j][i2 + dims_pad[2] * (i1 + dims_pad[1] * i0)]); */
        /*         printf("\n"); */
        /*     } */
        /*     printf("\n"); */
        /* } */
    }
    

    DEBUGPRINT("### Convolving");
    double delta = 2.0 * M_PI / dims[0];
    err = ambit_symmetric_odd_isotropic_circular_convolution_inplace(
        n_threads, n_k, k_abscissa, k_ordinate, dims[0], delta, x1, x2, x3);
    if (err) {
        printf("Error in ambit_symmetric_odd_isotropic_circular_convolution_inplace.\n");
        goto cleanup;
    }

    DEBUGPRINT("### Writing output");
    output_file_id = H5Fcreate(args.output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (output_file_id < 0) {
        printf("Error: Could not open \"%s\".\n", args.output_file);
        err = -1;
        goto cleanup;
    }
    output_dataspace_id = H5Screate_simple(4, dims, NULL);
    output_dataset_id   = H5Dcreate(output_file_id, "/simulation", H5T_NATIVE_DOUBLE, output_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int j = 0; j < 3; j++) {
        /* Define hyperslap in the file dataspace */
        offset[3] = j;
        count[3] = 1;
        err = H5Sselect_hyperslab(output_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        /* Write data to hyperslab */
        err = H5Dwrite(output_dataset_id, H5T_NATIVE_DOUBLE, 
                       memspace, output_dataspace_id,
                       H5P_DEFAULT, x[j]);
        if (err < 0) {
            printf("Error: Could not write hyperslab.\n");
            err = -1;
            goto cleanup;
        }
    }
    
  cleanup:
    if (memspace > 0) H5Sclose(memspace);
    if (output_dataspace_id > 0) H5Sclose(output_dataspace_id);
    if (output_dataset_id > 0) H5Dclose(output_dataset_id);
    if (output_file_id > 0) H5Fclose(output_file_id);
    if (levy_basis_dataspace_id > 0) H5Sclose(levy_basis_dataspace_id);
    if (levy_basis_dataset_id > 0) H5Dclose(levy_basis_dataset_id);
    if (levy_basis_file_id > 0) H5Fclose(levy_basis_file_id);
    if (kernel_file_id > 0) H5Fclose(kernel_file_id);
    free(tmp);
    free(k_abscissa);
    free(k_ordinate);
    free(x1);
    free(x2);
    free(x3);    
    return err;
}
