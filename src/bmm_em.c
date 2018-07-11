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

bmm_em_result em(Dataset* ds, int K, int max_iter, int verbose) {
  
  /** Allocate space **/
    
  //double** protos = sample_prototypes(ds, K);
  double** protos = sample_prototypes_hypercube(ds, K);
  double* pis = sample_pis(K);
  double** z = alloc_z(ds->N, K, ds->D);
  
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
    }
    
    
    // M-Step /////////////
    p_k(pis, z, K, ds->N);
    
    for (int k = 0; k < K; k++) {
      proto_k(ds, z, protos[k], k);
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
    
    // normalize by dividing by rowsums
    for (int k = 0; k < K; k++) {
      z[n][k] -= rowsum;
      z[n][k] = exp(z[n][k]);
      
      ll += z[n][k] * tmp[k];
  
    }
      
    free(tmp);
  }
    
  return ll;
  
}

void p_k(double* pis, double** z, int K, int N) {
  
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

double loglik(Dataset* ds, double** z, double* pis, double** protos, int K) {
  
  double ll = 0;
  
  for (int n = 0; n < ds->N; n++) {
    
    for (int k = 0; k < K; k++) {
      
      ll += z[n][k] * (log(pis[k]) + log_p_xn_k(ds, n, protos[k]));
      //Rprintf("z[n][k]= %f\n", z[n][k]);
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
