#pragma once
#ifndef ESTIMATEVARCOMP_H
#define ESTIMATEVARCOMP_H
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/unsupported/Eigen/MatrixFunctions"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

using namespace Eigen;

int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);
double cal_median(std::vector<double> col);
void cal_standard_deviation(std::vector<double>& ret, Eigen::MatrixXd& mat);
void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat, Eigen::MatrixXd& colmat);

Eigen::MatrixXd cor_mat(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, std::vector<double> vc_1, std::vector<double> vc_2);
Eigen::MatrixXd read_mat(std::ifstream& input_file, int row, int col);
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

#endif // ESTIMATEVARCOMP_H