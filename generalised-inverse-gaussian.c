#include "generalised-inverse-gaussian.h"

double gig_psi      (double lambda, double omega, double x) {
  return -lambda * (sinh(x) - x) - sqrt(lambda * lambda + omega * omega) * (cosh(x) - 1.0);
}

double gig_dpsi     (double lambda, double omega, double x) {
  return -lambda * (cosh(x) - 1.0) - sqrt(lambda * lambda + omega * omega) * sinh(x);
}

void   gig_psi_dpsi (double lambda, double omega, double x, double *y, double *dy) {
  double sq  = sqrt(lambda * lambda + omega * omega);
  double cm1 = cosh(x) - 1.0;
  double s   = sinh(x);
  *y  = -lambda * (s - x) - sq * cm1;
  *dy = -lambda * cm1 - sq * s;
} 

double gamma_psi      (double lambda, double x) {
  return -lambda * (expm1(x) - x);
}

double gamma_dpsi     (double lambda, double x) {
  return -lambda * expm1(x);
}

void   gamma_psi_dpsi (double lambda, double x, double *y, double *dy) {
  double em1 = expm1(x);
  *y  = -lambda * (em1 - x);
  *dy = -lambda * em1;
}

double gig_newton_raphson   (double lambda, double omega, double y0, double x0, 
                             double epsilon, int max_iter) {
  double x  = 0.0;
  double y  = 0.0;
  double dy = 0.0;
  
  for (int i = 0; i < max_iter; i++) {
    gig_psi_dpsi(lambda, omega, x0, &y, &dy);
    x = x0 - (y - y0) / dy;
    if (fabs(y - y0) < epsilon) break;
    x0 = x;
  }
  return x;
}

double gamma_newton_raphson (double lambda, double y0, double x0, double epsilon, int max_iter) {
  double x  = 0.0;
  double y  = 0.0;
  double dy = 0.0;
  
  for (int i = 0; i < max_iter; i++) {
    gamma_psi_dpsi(lambda, x0, &y, &dy);
    x = x0 - (y - y0) / dy;
    if (fabs(y - y0) < epsilon) break;
    x0 = x;
  }
  return x;
}

int gen_gig (int n_threads, gsl_rng **rng, 
             double lambda, double omega, 
             ptrdiff_t n, double *x) {
  if (lambda < 0 || omega <= 0) return -1;
  
  double t = gig_newton_raphson(lambda, omega, -1.0,  1.0, 0.01, 20);
  double s = gig_newton_raphson(lambda, omega, -1.0, -1.0, 0.01, 20);
  assert (t > 0 && s < 0);
  
  double psi_t, dpsi_t, psi_s, dpsi_s;
  gig_psi_dpsi(lambda, omega, t, &psi_t, &dpsi_t);
  gig_psi_dpsi(lambda, omega, s, &psi_s, &dpsi_s);

  double exp_mode = (lambda + sqrt(lambda * lambda + omega * omega)) / omega;
 
  double p, q, r, t1, s1;
  p  =  1.0 / dpsi_s;
  r  = -1.0 / dpsi_t;
  s1 = s - p * psi_s;
  t1 = t + r * psi_t;
  q  = t1 - s1;
  assert(p > 0 && q > 0 && r > 0 && t1 > 0 && s1 < 0);
  
  double   u, v, w, chi, z;
  int      thread_id;
  gsl_rng *rng0;
#pragma omp parallel private(thread_id, rng0, u, v, w, chi, z) shared(x) num_threads(n_threads)
  {
    thread_id = omp_get_thread_num();
    rng0 = rng[thread_id];
      
#pragma omp for 
    for (ptrdiff_t i = 0; i < n; i++) {
      /* printf(" %i, %ti\n", thread_id, i); */
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
          chi = 1.0;
        else if (z > t1)
          chi = exp(psi_t + dpsi_t * (z - t));
        else 
          chi = exp(psi_s + dpsi_s * (z - s));
        if (w * chi <= exp(gig_psi(lambda, omega, z))) break;
      }
      x[i] = exp(z) * exp_mode;
    }
  }    
  return 0;
}
  
int gen_gamma (int n_threads, gsl_rng **rng, 
               double lambda, 
               ptrdiff_t n, double *x) {
  if (lambda <= 0) return -1;
  
  double t = gamma_newton_raphson(lambda, -1.0,  1.0, 0.01, 20);
  double s = gamma_newton_raphson(lambda, -1.0, -1.0, 0.01, 20);
  assert (t > 0 && s < 0);
  
  double psi_t, dpsi_t, psi_s, dpsi_s;
  gamma_psi_dpsi(lambda, t, &psi_t, &dpsi_t);
  gamma_psi_dpsi(lambda, s, &psi_s, &dpsi_s);
  
  double exp_mode = lambda;

  double p, q, r, t1, s1;
  p  =  1.0 / dpsi_s;
  r  = -1.0 / dpsi_t;
  s1 = s - p * psi_s;
  t1 = t + r * psi_t;
  q  = t1 - s1;
  assert(p > 0 && q > 0 && r > 0 && t1 > 0 && s1 < 0);
  
  double   u, v, w, chi, z;
  int      thread_id;
  gsl_rng *rng0;
#pragma omp parallel private(thread_id, rng0, u, v, w, chi, z) shared(x) num_threads(n_threads)
  {
    thread_id = omp_get_thread_num();
    rng0 = rng[thread_id];
    
#pragma omp for 
    for (ptrdiff_t i = 0; i < n; i++) {
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
          chi = 1.0;
        else if (z > t1)
          chi = exp(psi_t + dpsi_t * (z - t));
        else 
          chi = exp(psi_s + dpsi_s * (z - s));
        if (w * chi <= exp(gamma_psi(lambda, z))) break;
      }
      x[i] = exp(z) * exp_mode;
    }
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
int generate_generalised_inverse_gaussian (int n_threads, gsl_rng **rng,
                                           double lambda, double chi, double psi, 
                                           ptrdiff_t n, double *x) {
  int err  = 0;
  double c = 0.0;
  /* Gamma */
  if (psi > 0.0 && chi == 0.0 && lambda > 0.0) {
    c = 2.0 / psi;
    err = gen_gamma(n_threads, rng, lambda, n, x);
    if (!err)
#pragma omp parallel for num_threads(n_threads)
      for (ptrdiff_t i = 0; i < n; i++)
        x[i] = c * x[i]; 
  }
  
  /* Inverse Gamma */
  else if (psi == 0.0 && chi > 0.0 && lambda < 0.0) {
    c = 2.0 / chi;
    err = gen_gamma(n_threads, rng, -lambda, n, x);
    if (!err)
#pragma omp parallel for num_threads(n_threads)
      for (ptrdiff_t i = 0; i < n; i++)
        x[i] = c / x[i]; 
  }
  
  /* GIG */
  else if (chi > 0.0 && psi > 0.0 && lambda >= 0.0) {
    c = sqrt(chi / psi);
    err = gen_gig(n_threads, rng, lambda, sqrt(chi * psi), n, x);
    if (!err)
#pragma omp parallel for num_threads(n_threads)
      for (ptrdiff_t i = 0; i < n; i++)
        x[i] = c * x[i];
  }
  
  /* Inverse GIG */
  else if (psi > 0.0 && chi > 0.0 && lambda < 0.0) {
    c = sqrt(chi / psi);
    err = gen_gig(n_threads, rng, -lambda, sqrt(chi * psi), n, x);
    if (!err) 
#pragma omp parallel for num_threads(n_threads)
      for (ptrdiff_t i = 0; i < n; i++)
        x[i] = c / x[i];
  }
  
  /* Invalid parameters */
  else
    err = -1;
  
  return err;
}
