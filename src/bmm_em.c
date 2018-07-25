#include "bmm_em.h"
#include "float.h"
#include "Rinternals.h"
#include "R.h"
#include <omp.h>


void free_bmm_em_result(bmm_em_result* x) {
  
  free(x->pis);
  for (int k = 0; k < x->K; k++) {
    free(x->protos[k]);
  }
  free(x->protos);
  free(x->cluster);
  
}

bmm_em_result em(Dataset* ds, int K, int max_iter, int verbose, int hbbmm) {
  
  /** Allocate space **/
    
  //double** protos = sample_prototypes(ds, K);
  double** protos = sample_prototypes_hypercube(ds, K);
  double* pis = sample_pis(K);
  double** z = alloc_z(ds->N, K, ds->D);
  
  // calculate alpha and beta estimates
  double alpha, beta;
  empirical_bayes(ds, &alpha, &beta);
  Rprintf("alph=%f | beta=%f\n", alpha, beta);

  double thresh = 1e-6;
  int converged = 0;
  double prev = 0;
  double ll = -DBL_MAX;
  int iter = 0;
  
  while (iter < max_iter && !converged) {
    prev = ll;
    ll = 0;
    
    // Expectation & log likelihood
    ll = log_z_nk(ds, z, pis, protos, K);
    
    
    if (verbose) {
      Rprintf(" %4d | %15.4f\n", iter, ll);
    }
    
    // Check converged
    if (ll - prev < thresh) {
      converged = 1;
      if (verbose) {
        Rprintf("-- Converged --\n");
      }
      break;
    }
    
    
    // M-Step /////////////
    p_k(pis, z, K, ds->N);
    
    // Pass alpha and beta here
    for (int k = 0; k < K; k++) {
      
      if (hbbmm) {
        proto_k_hbbmm(ds, z, protos[k], k, alpha, beta);
      } else {
        proto_k(ds, z, protos[k], k);
      }
    }
    // End M-Step /////////
    
    iter++;
  }
  
  // get cluster
  int * cluster = (int*) calloc(ds->N, sizeof(int));
  
  for (int n = 0; n < ds->N; n++) {
    
    double max = 0;
    
    for (int k = 0; k < K; k++) {
      if (z[n][k] > max) {
        max = z[n][k];
        cluster[n] = k;
      }
    }
  }
  
  // free everything
  for (int n = 0; n < ds->N; n++) {
    free(z[n]);
  }
  free(z);
  
  bmm_em_result result = {
    .protos = protos,
    .pis = pis,
    .cluster = cluster,
    .D = ds->D,
    .K = K,
    .ll = ll
  };
  
  return result;
  
}

double clip(double x) {
  double lo = 0.000000001;
  double hi = 0.999999999;
  
  if (x < lo) {
    x = lo;
  } else if (x > hi) {
    x = hi;
  }
  return x;
}

// likelihood of xn given cluster k
double log_p_xn_k(Dataset* ds, int row, double* proto) {
  
  double ll = 0;
  double mu;
  
  for (int i = 0; i < ds->D; i++) {
    
    mu = clip(proto[i]);
    //mu = proto[i];
    
    ll += at(ds, row, i) ? log(mu) : log(1 - mu);
  }
  return ll;
  
}

double log_z_nk(Dataset* ds, double** z, double* pis, double** protos, int K) {
  
  double ll = 0; 
  double * tmp;
  
  #pragma omp parallel for private(tmp) reduction (+:ll)
  for (int n = 0; n < ds->N; n++) {
    
    tmp = (double*) malloc(K * sizeof(double));
    double max = -DBL_MAX;
    
    for (int k = 0; k < K; k++) {
      
      tmp[k] = log(pis[k]) + log_p_xn_k(ds, n, protos[k]);
      
      z[n][k] = tmp[k];
    
      // Keep track of the max for logsumexp trick  
      if (z[n][k] > max) {
        max = z[n][k];
      }
    }
    
    // logsumexp trick
    double rowsum = 0;
    for (int k = 0; k < K; k++) {
      rowsum += exp(z[n][k] - max);
    }
    
    rowsum = max + log(rowsum);
    
    //Rprintf("z[n][k]=");
    
    // normalize by dividing by rowsums
    for (int k = 0; k < K; k++) {
      z[n][k] -= rowsum;
      z[n][k] = exp(z[n][k]);
      
      //Rprintf(" %f", z[n][k]);
      
      ll += z[n][k] * tmp[k];
  
    }
    //Rprintf("\n");
      
    free(tmp);
  }
    
  return ll;
  
}

void p_k(double* pis, double** z, int K, int N) {
  
  #pragma omp parallel for shared(pis, z)
  for (int k = 0; k < K; k++) {
    pis[k] = 0; // reset to zero
    
    for (int n = 0; n < N; n++) {
      
      pis[k] += z[n][k];
      
    }
    
    // make sure pis are not zero or one
    pis[k] /= N;
    
  }
  
}

void proto_k(Dataset* ds, double** z, double* proto, int k) {
  
  #pragma omp parallel for shared(proto, z)
  for (int i = 0; i < ds->D; i++) {
    
    double num = 0;
    double den = 0;
    
    for (int n = 0; n < ds->N; n++) {
      
      num += z[n][k] * at(ds, n, i);
      den += z[n][k];
      
    }
    
    proto[i] = clip(num / den);
    
  }
}

void proto_k_hbbmm(Dataset* ds, double** z, double* proto, int k, double alpha, double beta) {
  
  #pragma omp parallel for shared(proto, z)
  for (int i = 0; i < ds->D; i++) {
    
    double num = 0;
    double den = 0;
    
    for (int n = 0; n < ds->N; n++) {
      
      num += z[n][k] * at(ds, n, i);
      den += z[n][k];
      
    }
    
    //proto[i] = clip( (num  / den);
    proto[i] = (num + alpha - 1.0) / (den + alpha + beta - 2.0);
    
  }
}


double loglik(Dataset* ds, double** z, double* pis, double** protos, int K) {
  
  double ll = 0;
  
  for (int n = 0; n < ds->N; n++) {
    
    for (int k = 0; k < K; k++) {
      
      ll += z[n][k] * (log(pis[k]) + log_p_xn_k(ds, n, protos[k]));
      ////Rprintf("z[n][k]= %f\n", z[n][k]);
    }
  }
  
  return ll;
  
}

double* sample_pis(int K) {
  
  double* pis = (double*) calloc(K, sizeof(double));
  
  for (int k = 0; k < K; k++) {
    pis[k] = 1.0/K;
  }
  return pis;
}


double** sample_prototypes_hypercube(Dataset* ds, int K)  {
  
  double** protos = (double**) calloc(K, sizeof(double*));
  for (int k = 0; k < K; k++) {
    protos[k] = (double*) calloc(ds->D, sizeof(double));
  }
  
  // randomly sample a row from x 
  for (int k = 0; k < K; k++) {
    
    // loop over bits
    for (int i = 0; i < ds->D; i++) {
      GetRNGstate();
      double rand = (unif_rand() - 0.50) * 1e-2;
      PutRNGstate();
      
      protos[k][i] = 0.50 + rand;
      
    }
  }
  return protos;
}

double** sample_prototypes(Dataset* ds, int K)  {
  
  double** protos = (double**) calloc(K, sizeof(double*));
  for (int k = 0; k < K; k++) {
    protos[k] = (double*) calloc(ds->D, sizeof(double));
  }
  
  // randomly sample a row from x 
  for (int k = 0; k < K; k++) {
    
    GetRNGstate();
    int row = floor(unif_rand() * ds->N);
    PutRNGstate();
    
    // loop over bits
    for (int i = 0; i < ds->D; i++) {
      GetRNGstate();
      double rand = unif_rand();
      PutRNGstate();
      
      protos[k][i] = 0.25*at(ds, row, i) + 0.75*rand;
      
    }
  }
  return protos;
}

double** alloc_z(int N, int K, int D) {
  
  double** z = (double**) calloc(N, sizeof(double*));
  for (int n = 0; n < N; n++) {
    z[n] = (double*) calloc(K, sizeof(double));
    
    for (int k = 0; k < K; k++) {
      z[n][k] = 1.0/K;
    }
  }
  
  return z;
}


SEXP convert_bmm_em_result(Dataset* ds, bmm_em_result * res, int* prtCnt) {
  
  /** Return the following:
   * prototypes as matrix
   * pis as numeric vector
   * cluster as integer vector
   * K, D, ll?
   */ 
  
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  (*prtCnt)++;
  
  SEXP protos = PROTECT(Rf_allocMatrix(REALSXP, ds->D, res->K));
  (*prtCnt)++;
  
  SEXP pis = PROTECT(Rf_allocVector(REALSXP, res->K));
  (*prtCnt)++;
  
  SEXP cluster = PROTECT(Rf_allocVector(INTSXP, ds->N));
  (*prtCnt)++;
  
  // TODO: change iteration to avoid tranposing on back end
  // copy data to R vectors
  for (int k = 0; k < res->K; k++) {
    memcpy(&REAL(protos)[k * ds->D], res->protos[k], ds->D * sizeof(double));
  }
  
  memcpy(REAL(pis), res->pis, res->K * sizeof(double));
  memcpy(INTEGER(cluster), res->cluster, ds->N * sizeof(int));
  
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  (*prtCnt)++;
  
  SET_VECTOR_ELT(out, 0, protos);
  SET_VECTOR_ELT(out, 1, pis);
  SET_VECTOR_ELT(out, 2, cluster);
  
  SET_STRING_ELT(names, 0, Rf_mkChar("prototypes"));
  SET_STRING_ELT(names, 1, Rf_mkChar("pis"));
  SET_STRING_ELT(names, 2, Rf_mkChar("cluster"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  
  return out;
  
}

znk_result predict_log_z_nk(Dataset *ds, double* pis, double** protos, int K) {
  
  double** z = alloc_z(ds->N, K, ds->D);
  
  double ll = 0;
  
  #pragma omp parallel for
  for (int n = 0; n < ds->N; n++) {
    
    for (int k = 0; k < K; k++) {
      
      z[n][k] = log_p_xn_k(ds, n, protos[k]);
      
    }
    
  }

  znk_result res = {.z = z, .ll = ll, .K = K};
  
  return (res);
  
};

void free_znk_result(znk_result* x) {
  for (int k = 0; k < x->K; k++) {
    free(x->z[k]);
  }
  free(x->z);
};


SEXP convert_znk_result(Dataset* ds, znk_result * res, int* prtCnt) {
  
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  (*prtCnt)++;
  
  SEXP znk = PROTECT(Rf_allocMatrix(REALSXP, ds->N, res->K));
  (*prtCnt)++;
  
  SEXP ll = PROTECT(Rf_allocVector(REALSXP, 1));
  (*prtCnt)++;
  
  // TODO: change iteration to avoid tranposing on back end
  // copy data to R vectors
  for (int n = 0; n < ds->N; n++) {
    for (int k = 0; k < res->K; k++) {
    REAL(znk)[k * ds->N + n] = res->z[n][k];
    }
  }
  
  REAL(ll)[0] = res->ll;
  
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  (*prtCnt)++;
  
  SET_VECTOR_ELT(out, 0, znk);
  SET_VECTOR_ELT(out, 1, ll);
  
  SET_STRING_ELT(names, 0, Rf_mkChar("z"));
  SET_STRING_ELT(names, 1, Rf_mkChar("ll"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  
  return out;
  
}

 
// functions for calculating beta params

double beta_hat(double N, double C, double a0) {
  double nca = (N - C)  * a0;
  double bhat = ( nca + C + sqrt( (nca + C)*(nca + C) - 4.0*C*nca) ) / (2.0*C);
  return bhat;
}

double alpha_hat(double N, double C, double bhat) {
  double cbn = (C*bhat + N);
  double ahat = (cbn + sqrt( cbn*cbn + 4.0*C*bhat*(C-N)) ) / (2.0*(N-C));
  return ahat;
}

double beta_nought(double N, double C) {
  return (N + sqrt(N*N - 4.0*C*(N-C)) )/(2.0*C);
}

double alpha_nought(double N, double C, double b0) {
  double bcn = (b0*C + N);
  double res = ( (b0*C - C + N) + sqrt( (bcn*bcn) + 4.0*C*(C-N)*b0) ) / (2.0*(N-C));
  return res;
}

void empirical_bayes(Dataset* ds, double* alpha, double* beta) {
  
  double C = 0.0;
  double N = (double) ds->N * ds->D;
  
  for (int n = 0; n < ds->N; n++) {
    for (int d = 0; d < ds->D; d++) {
      C += at(ds, n, d);
    }
  }
  
  double b0 = beta_nought(N, C);
  double a0 = alpha_nought(N, C, b0);
  double bhat = beta_hat(N, C, a0);
  double ahat = alpha_hat(N, C, bhat);

  *alpha = ahat;
  *beta = bhat;
}



