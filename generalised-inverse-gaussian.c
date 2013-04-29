#include "ambit-stochastics.h"

//
// Utility functions
//

inline R
gig_psi(R lambda, R omega, R x) {
    return -lambda * (sinh(x) - x) - sqrt(lambda * lambda + omega * omega) * (cosh(x) - 1.0);
}

inline R 
gig_dpsi(R lambda, R omega, R x) {
    return -lambda * (cosh(x) - 1.0) - sqrt(lambda * lambda + omega * omega) * sinh(x);
}

inline void 
gig_psi_dpsi(R lambda, R omega, R x, R *y, R *dy) {
    R sq  = sqrt(lambda * lambda + omega * omega);
    R cm1 = cosh(x) - 1.0;
    R s   = sinh(x);
    *y  = -lambda * (s - x) - sq * cm1;
    *dy = -lambda * cm1 - sq * s;
} 

inline R 
gamma_psi(R lambda, R x) {
    return -lambda * (expm1(x) - x);
}

inline R 
gamma_dpsi(R lambda, R x) {
    return -lambda * expm1(x);
}

inline void 
gamma_psi_dpsi(R lambda, R x, R *y, R *dy) {
    R em1 = expm1(x);
    *y  = -lambda * (em1 - x);
    *dy = -lambda * em1;
}

R 
gig_newton_raphson(R lambda, R omega, R y0, R x0, 
                   R epsilon, int max_iter) {
    R x  = 0.0;
    R y  = 0.0;
    R dy = 0.0;
  
    for (int i = 0; i < max_iter; ++i) {
        gig_psi_dpsi(lambda, omega, x0, &y, &dy);
        x = x0 - (y - y0) / dy;
        if (fabs(y - y0) < epsilon) break;
        x0 = x;
    }
    return x;
}

R 
gamma_newton_raphson(R lambda, R y0, R x0, R epsilon, int max_iter) {
    R x  = 0.0;
    R y  = 0.0;
    R dy = 0.0;
  
    for (int i = 0; i < max_iter; ++i) {
        gamma_psi_dpsi(lambda, x0, &y, &dy);
        x = x0 - (y - y0) / dy;
        if (fabs(y - y0) < epsilon) break;
        x0 = x;
    }
    return x;
}

//
// Generators
//

int
gen_gig (gsl_rng **rng, struct ambit_dense_array *a, R lambda, R chi, R psi) {
    bool inverse;
    R c;

    // GIG
    if (chi > 0.0 && psi > 0.0 && lambda >= 0.0) {
        c = sqrt(chi / psi);
        inverse = false;
    }
    // Inverse GIG
    else if (psi > 0.0 && chi > 0.0 && lambda < 0.0) {
        c = sqrt(chi / psi);
        lambda = -lambda;
        inverse = true;
    }
    else
        return -1;
    
    R omega = sqrt(chi * psi);

    R t = gig_newton_raphson(lambda, omega, -1.0,  1.0, 0.01, 20);
    R s = gig_newton_raphson(lambda, omega, -1.0, -1.0, 0.01, 20);
    assert (t > 0 && s < 0);
  
    R psi_t, dpsi_t, psi_s, dpsi_s;
    gig_psi_dpsi(lambda, omega, t, &psi_t, &dpsi_t);
    gig_psi_dpsi(lambda, omega, s, &psi_s, &dpsi_s);

    R c_exp_mode         = c * (lambda + sqrt(lambda * lambda + omega * omega)) / omega;
    R c_exp_mode_inverse = c * omega / (lambda + sqrt(lambda * lambda + omega * omega));

    R p, q, r, t1, s1;
    p  =  1.0 / dpsi_s;
    r  = -1.0 / dpsi_t;
    s1 = s - p * psi_s;
    t1 = t + r * psi_t;
    q  = t1 - s1;
    assert(p > 0 && q > 0 && r > 0 && t1 > 0 && s1 < 0);
  
    R  u, v, w, chi0, z;
    gsl_rng *rng0;
    size_t *sub0 = NULL;
    size_t lin_idx0;
    size_t inner_n = a->n / a->dim[0];
#pragma omp parallel                                \
    private(rng0, sub0, lin_idx0, u, v, w, chi0, z) \
    shared(a)
    {
        rng0 = rng[omp_get_thread_num()];
        sub0 = malloc(a->rank * sizeof(size_t));
#pragma omp for 
        for (size_t i0 = 0; i0 < a->dim[0]; ++i0) {
            sub0[0] = i0;
            for (int j = 1; j < a->rank; ++j)
                sub0[j] = 0;
            for (size_t i = 0; i < inner_n; ++i) {
                lin_idx0 = ambit_dense_array_linear_index(a, sub0);
                while (true) {
                    u = gsl_rng_uniform(rng0);
                    v = gsl_rng_uniform(rng0);
                    w = gsl_rng_uniform(rng0);
                    if (v == 0.0) continue;
                    if (u < q / (p + q + r)) 
                        z = s1 + q * v;
                    else if (u < (q + r) / (p + q + r))
                        z = t1 - r * log(v);
                    else 
                        z = s1 + p * log(v);
                    if (s1 <= z && z <= t1)
                        chi0 = 1.0;
                    else if (z > t1)
                        chi0 = exp(psi_t + dpsi_t * (z - t));
                    else 
                        chi0 = exp(psi_s + dpsi_s * (z - s));
                    if (w * chi0 <= exp(gig_psi(lambda, omega, z))) break;
                }
                a->data[lin_idx0] = inverse ? exp(-z) * c_exp_mode_inverse : exp(z) * c_exp_mode;
                increment_subscript(sub0, a->rank, a->dim);
            }
        }
        free(sub0);
    }    
    return 0;
}
  
int 
gen_gamma (gsl_rng **rng, struct ambit_dense_array *a, R lambda, R c) {
     
    if (c <= 0 || lambda == 0)
        return -1;
    
    bool inverse = lambda < 0;
    lambda = fabs(lambda);

    R t = gamma_newton_raphson(lambda, -1.0,  1.0, 0.01, 20);
    R s = gamma_newton_raphson(lambda, -1.0, -1.0, 0.01, 20);
    assert (t > 0 && s < 0);
  
    R psi_t, dpsi_t, psi_s, dpsi_s;
    gamma_psi_dpsi(lambda, t, &psi_t, &dpsi_t);
    gamma_psi_dpsi(lambda, s, &psi_s, &dpsi_s);
  
    R c_exp_mode = c * lambda;
    R c_exp_mode_inverse = c / lambda;

    R p, q, r, t1, s1;
    p  =  1.0 / dpsi_s;
    r  = -1.0 / dpsi_t;
    s1 = s - p * psi_s;
    t1 = t + r * psi_t;
    q  = t1 - s1;
    assert(p > 0 && q > 0 && r > 0 && t1 > 0 && s1 < 0);
  
    R   u, v, w, chi0, z;
    gsl_rng *rng0;
    size_t *sub0 = NULL;
    size_t lin_idx0;
    size_t inner_n = a->n / a->dim[0];
#pragma omp parallel                                           \
    private(rng0, sub0, lin_idx0, u, v, w, chi0, z) \
    shared(a) 
    {
        rng0 = rng[omp_get_thread_num()];
        sub0 = malloc(a->rank * sizeof(size_t));
#pragma omp for 
        for (size_t i0 = 0; i0 < a->dim[0]; ++i0) {
            sub0[0] = i0;
            for (int j = 1; j < a->rank; ++j)
                sub0[j] = 0;
            for (size_t i = 0; i < inner_n; ++i) {
                lin_idx0 = ambit_dense_array_linear_index(a, sub0);
                while (true) {
                    u = gsl_rng_uniform(rng0);
                    v = gsl_rng_uniform(rng0);
                    w = gsl_rng_uniform(rng0);
                    if (v == 0.0) continue;
                    if (u < q / (p + q + r)) 
                        z = s1 + q * v;
                    else if (u < (q + r) / (p + q + r))
                        z = t1 - r * log(v);
                    else 
                        z = s1 + p * log(v);
                    if (s1 <= z && z <= t1)
                        chi0 = 1.0;
                    else if (z > t1)
                        chi0 = exp(psi_t + dpsi_t * (z - t));
                    else 
                        chi0 = exp(psi_s + dpsi_s * (z - s));
                    if (w * chi0 <= exp(gamma_psi(lambda, z))) break;
                }
                a->data[lin_idx0] = inverse ? exp(-z) * c_exp_mode_inverse : exp(z) * c_exp_mode;
                increment_subscript(sub0, a->rank, a->dim);
            }
        }
        free(sub0);
    }    
    return 0;
}

/** Generate generalised inverse Gaussian random variates
 *
 * n_threads: number of OpenMP threads to use
 * rng: random number generators, array of length n_threads
 * lambda, chi, psi: GIG-parameters corresponding to a density proportional to 
 *                   (psi/chi)^(lambda/2) x^(lambda-1) exp(-1/2 (chi/x + psi x))
 * n: number of variates to generate
 * x: array of length n to store the result 
 */
int 
generate_generalised_inverse_gaussian (gsl_rng **rng, struct ambit_dense_array *a, R lambda, R chi, R psi) {
    int err  = 0;
    R inf = atof("inf");
    /* Dirac delta.
     * The limit of GIG(lambda, c^2 z, z) as z goes to infinity is delta_c.
     * We represent this limit as GIG(lambda, c^2, inf).
     */
    if (chi >= 0 && psi == inf)
        ambit_dense_array_set(a, sqrt(chi));
    
    // Gamma and Inverse Gamma
    else if ((psi > 0.0 && chi == 0.0 && lambda > 0.0) || (psi == 0.0 && chi > 0.0 && lambda < 0.0)) 
        err = gen_gamma(rng, a, lambda, 2.0 / fmax(psi, chi));
  
    // GIG and Inverse GIG
    else if (chi > 0.0 && psi > 0.0) 
        err = gen_gig(rng, a, lambda, chi, psi); 

    // Invalid parameters
    else
        err = -1;
  
    return err;
}
