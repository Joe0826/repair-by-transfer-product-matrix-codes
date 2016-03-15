#ifndef MBR_MSR_H
#define MBR_MSR_H
#include "product_matrix.h"

GeneratorDAG* make_MBR_decoder_dag(int* psi, int* M, int* nodes, int k, int d, int w);

GeneratorDAG* make_MSR_decoder_dag(int* psi, int* lambdas, int* M, int* nodes, int k, int d, int w);

GeneratorDAG* make_MSR_2k2_decoder_dag(GeneratorDAG* src, int* psi, int* lambdas, int* M, int* nodes, int k, int w);

GeneratorPlan* make_MSR_2k2_decoder_plan(int* psi, int* lambdas, int* M, int* nodes, int k, int w);

/*
  Generate matrix M
*/
int* make_m_MBR(int k, int d);

int* make_psi_MBR(int n, int k, int d, int w);

void make_mu_MBR_inplace(int* mu, const int f, const int k, const int d, const int w);

/*
  Generate mu_f (the helper-vector for node f failing).
*/
int* make_mu_MBR(int f, int k, int d, int w);

/*
Generate d-by-(d - k + 1) matrix
*/
int* make_m_MSR(int k, int d);

int* make_psi_MSR(int n, int k, int d, int w);

void make_psi_MSR_ext(int* psi, int* lambdas, int n, int k, int d, int w);

void make_mu_MSR_inplace(int* mu, const int f, const int k, const int d, const int w);

/*
  Generate mu_f (the helper-vector for node f failing).
*/
int* make_mu_MSR(int f, int k, int d, int w);

#endif
