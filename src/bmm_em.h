#ifndef BMM_EM_H
#define BMM_EM_H

#include "dataset.h"
#include "R.h"
#include "Rinternals.h"

typedef struct {
  double ** protos;
  double* pis;
  int * cluster;
  int K;
  int D;
  double ll;
} bmm_em_result;

typedef struct {
  double ** z;
  double ll;
  int K;
} znk_result;

void free_bmm_em_result(bmm_em_result* x);

bmm_em_result em(Dataset* ds, int K, int max_iter, int verbose, int hbbmm);

double clip(double x);

double log_p_xn_k(Dataset* ds, int row, double* proto);

double log_z_nk(Dataset* ds, double** znk, double* pis, double** protos, int K);

void p_k(double* pis, double** z, int K, int N);

void proto_k(Dataset* ds, double** z, double* proto, int k);

void proto_k_hbbmm(Dataset* ds, double** z, double* proto, int k, double alpha, double beta);

double loglik(Dataset* ds, double** z, double* pis, double** protos, int K);

double* sample_pis(int K);

double** sample_prototypes(Dataset* ds, int K);

double** sample_prototypes_hypercube(Dataset* ds, int K);

double** alloc_z(int N, int K, int D);

znk_result predict_log_z_nk(Dataset *ds, double* pis, double** protos, int K);

void free_znk_result(znk_result* x);

SEXP convert_bmm_em_result(Dataset* ds, bmm_em_result * res, int * prtCnt);

SEXP convert_znk_result(Dataset* ds, znk_result* res, int * prtCnt);

double beta_hat(double N, double C, double a0);

double alpha_hat(double N, double C, double bhat);

double beta_nought(double N, double C);

double alpha_nought(double N, double C, double b0);

void empirical_bayes(Dataset* ds, double* alpha, double* beta);

#endif /* BMM_EM_H */

