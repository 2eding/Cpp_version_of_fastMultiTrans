#pragma once

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/unsupported/Eigen/MatrixFunctions"
#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>


#define _FILE_OFFSET_BITS  64 // large file support

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mkl.h>

#define MAXNSEC 1000 // maximum number of sections for scaling factor
#define EIGEN_EPS 1E-10 // epsilon to add to the matrix diagonal
#define DBL_EPS_HARD 1E-20 // hard epsilon
#define DBL_EPS_SOFT 1E-10 // soft epsilon
#define P_THRES .05 // starting point of scaling factor computation
#define BONFFACTOR 10 // conservative estimate of bonferroni factor
#define DECREASE_RATIO .1 // control between-sections of scaling factor part
#define DOSE_EFF .5 // Armitage test
#define ONEWRITE 1000
#define DBL_CHUNK 1000
#define EPS 1E-100

#define min(a,b) (((a)<(b)) ? (a) : (b))
#define max(a,b) (((a)>(b)) ? (a) : (b))

using namespace Eigen;

int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);
double cal_median(std::vector<double> col);
void cal_standard_deviation(std::vector<double>& ret, Eigen::MatrixXd& mat);
void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat, Eigen::MatrixXd& colmat);

Eigen::MatrixXd cor_mat(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, std::vector<double> vc_1, std::vector<double> vc_2);
Eigen::MatrixXd read_mat(std::ifstream& input_file, int row, int col);
Eigen::MatrixXd corrband_mat(int windowSize, Eigen::MatrixXd cor_mat);
Eigen::MatrixXd estimateKinship(Eigen::MatrixXd snp);
void estimateVarComp(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, Eigen::MatrixXd phenotype, std::vector<double>& vc_1, std::vector<double>& vc_2);
struct lmm_fit {
    double hsq;
    VectorXd beta;
    double sigmasq;
    double loglik;
    double rss;
    double logdetXSX;
};

struct calcLL_args {
    VectorXd Kva;
    VectorXd y;
    MatrixXd X;
    bool reml;
    double logdetXpX;
};

struct eigenrot {
    VectorXd Kva;
    MatrixXd Kve;
    MatrixXd y;
    MatrixXd X;
};

// calc X'X
MatrixXd calc_XpX(const MatrixXd& X);

// eigen decomposition
//    returns eigenvalues and transposed eigenvectors
std::pair<VectorXd, MatrixXd> eigen_decomp(const MatrixXd& A);

// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(const MatrixXd& K, const MatrixXd& y, const MatrixXd& X);

// calculate log det X'X
double calc_logdetXpX(const MatrixXd& X);

// getMLsoln
// for fixed value of hsq, calculate MLEs of beta and sigmasq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = whether you'll be using REML (so need to calculate log det XSX)
struct lmm_fit getMLsoln(const double hsq, const VectorXd& Kva, const VectorXd& y,
    const MatrixXd& X, const bool reml);

// calcLL
// calculate log likelihood for fixed value of hsq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// logdetXpX = log det X'X; if NA, it's calculated
struct lmm_fit calcLL(const double hsq, const VectorXd& Kva, const VectorXd& y,
    const MatrixXd& X, const bool reml, const double logdetXpX);

// just the negative log likelihood, for the optimization
double negLL(const double x, struct calcLL_args* args);

// fitLMM
// Optimize log liklihood over hsq
//
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// check_boundary = if true, explicity check 0.0 and 1.0 boundaries
// logdetXpX = log det X'X; if NA, it's calculated
// tol   = tolerance for convergence
struct lmm_fit fitLMM(const VectorXd& Kva, const VectorXd& y, const MatrixXd& X,
    const bool reml, const bool check_boundary,
    const double logdetXpX, const double tol);

// 1-d optimization by Brent's method
double qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void*),
    void* info, double tol);

//slide.h
int compare_by_stat(const void* a, const void* b);
int compare_double(const void* s, const void* t);
int compare_double_rev(const void* s, const void* t);

// for analysis
void correct_pvalue(double*, int, double, double*, double*);
void per_marker_threshold(double*, int, double, double*, double*);

FILE* readfile(char* filename);
FILE* writefile(char* filename);

void slide_1prep(int windowsize, char* input_c, char* output);
void slide_2run(char* input, char* output);
void slide_3sort(char* input_max, char* output_sor);
void slide_4correct(char* input_sor, char* output);