#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

struct arguments;
int validate_arguments(struct arguments *p);


const char *argp_program_version     = "simulate-trawl-process v. 0.0";
const char *argp_program_bug_address = "<emil@math.au.dk>";

/* Program documentation. */
static char doc[] = "A program to simulate time series from trawl processes";

/* A description of the arguments we accept. */
static char args_doc[] = "n output";

/* The options we understand. */
static struct argp_option options[] = {
    {0, 0, 0, OPTION_DOC, "Precalculated ambit set:"},
    {"ambit-set",            'A', "file",   0, "Use the ambit set specified in 'file'", 1},
    {0, 0, 0, OPTION_DOC, "Parametrised cascade ambit set:"},   
    {"cascade",              'C', 0,        0, "Use the \"cascade\" ambit set",                 2 },
    {"decorelation-length",  'T', "T",      0, "Decorrelation length of the cascade ambit set", 2 },
    {"cascade-length",       'L', "L",      0, "Cascade length of the cascade ambit set",       2 },
    {"smoothness",           't', "t",      0, "Smoothness parameter of the cascade ambit set", 2 },
    {"resolution",           'r', "dt",     0, "Simulation step size",                          2 },
    {0, 0, 0, OPTION_DOC, "Levy seed:"},
    {"normal",               'N', 0,        0, "Use a normal seed", 3 },
    {"generalised-hyperbolic", 'H', 0,        0, "Use a generalised hyperbolic seed", 3 },
    {"lambda",               'l', "lambda", 0, "Parameter of GH seed", 3 },
    {"alpha",                'a', "alpha",  0, "Parameter of GH seed", 3 },
    {"beta",                 'b', "beta",   0, "Parameter of GH seed", 3 },
    {"mu",                   'm', "mu",     0, "Parameter of GH seed or mean of normal seed", 3 },
    {"delta",                'd', "delta",  0, "Parameter of GH seed or variance of normal seed", 3 },
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
    case 't': // theta
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
    case ARGP_KEY_ARG:
        switch (state->arg_num) {
        case 0:
            arguments->n = atol(arg);
            break;
        case 1:
            arguments->output_file = arg;
            break;
        default:
            /* Too many arguments. */
            argp_usage (state);
        }            
        break;
    case ARGP_KEY_END:
        if (state->arg_num < 2 || validate_arguments(arguments))
            argp_usage (state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };


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
        printf("The simulation length 'n' must be a positive integer.\n");
        return -1;
    }
    if (!p->output_file) {
        printf("Please specify an output file.\n");
        return -1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    struct arguments arguments;
    initialise_arguments(&arguments);
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    printf_arguments(&arguments);
    return 0;
}
