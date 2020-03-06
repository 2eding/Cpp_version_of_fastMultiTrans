#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "estimateVarComp.h"
#include <vector>
#include <list>
#define MAXBUFSIZE  ((int) 1e6)

using namespace Eigen;

int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);

//estimateVarComp input kinship, pheno, SNP, outputPath
int main() {
    //MatrixXd kin_mat = load_file<MatrixXd>(argv[1]); //Kinship matrix
    //MatrixXd pheno_mat = load_file<MatrixXd>(argv[2]); //Phenotype matrix
    //MatrixXd snp_mat = load_file<MatrixXd>(argv[3]); //SNP matrix
    //    
    /*std::string kin = ;
    std::string pheno = argv[2];
    std::string snp = argv[3];*/

    std::ifstream kinin("K.txt");
    std::ifstream phenoin("Y.txt");
    std::ifstream snpin("X.txt");
    std::ofstream out("VC2.txt");

    std::string read_buffer; //file preprocessing variable
    std::string token;
    std::stringstream stream;


    int kin_mat_row(0), kin_mat_col(0);
    int phe_mat_row(0), phe_mat_col(0);
    int snp_mat_row(0), snp_mat_col(0);
    
    kin_mat_row = count_matrix_row(kinin); kin_mat_col = count_matrix_col(kinin);
    phe_mat_row = count_matrix_row(phenoin); phe_mat_col = count_matrix_col(phenoin);
    snp_mat_row = count_matrix_row(snpin); snp_mat_col = count_matrix_col(snpin);

    MatrixXd kin_mat = MatrixXd(kin_mat_row, kin_mat_col); //Kinship matrix
    MatrixXd pheno_mat = MatrixXd(phe_mat_row, phe_mat_col); //Phenotype matrix
    MatrixXd snp_mat = MatrixXd(snp_mat_row, snp_mat_col); //SNP matrix

    for (int i = 0; i < kin_mat_row; i++) { //row
        std::getline(kinin, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < kin_mat_col; j++) { //col
            stream >> token;
            kin_mat(i, j) = std::stold(token);

            //std::cout.precision(16); 
            //std::cout << X_right_matrix(i,j) << " ";
        }
        stream.clear();
        //std::cout << "" << std::endl;
    }
    for (int i = 0; i < phe_mat_row; i++) { //row
        std::getline(phenoin, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < phe_mat_col; j++) { //col
            stream >> token;
            pheno_mat(i, j) = std::stold(token);
            //std::cout.precision(16); 
            //std::cout << X_right_matrix(i,j) << " ";
        }
        stream.clear();
        //std::cout << "" << std::endl;
    }
    for (int i = 0; i < snp_mat_row; i++) { //row
        std::getline(snpin, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < snp_mat_col; j++) { //col
            stream >> token;
            snp_mat(i, j) = std::stold(token);
            //std::cout.precision(16); 
            //std::cout << X_right_matrix(i,j) << " ";
        }
        stream.clear();
        //std::cout << "" << std::endl;
    }
 

    if (true) {
        MatrixXd oriKin = kin_mat;
        struct eigenrot e;
        struct lmm_fit vc;
        for (int i = 1; i <= pheno_mat.rows(); i++) {
            kin_mat = oriKin;
            std::cout << i << std::endl;
            
            e = eigen_rotation(kin_mat, pheno_mat.row(i), snp_mat);
            vc = fitLMM(e.Kva, e.y, e.X, true, true, NULL, 1e-6);

            // vc.hsq * vc.sigmasq = Vg
            // (1 - vc.hsq) * vc.sigmasq = Ve
            out << vc.hsq * vc.sigmasq << "\t" << (1 - vc.hsq) * vc.sigmasq << "\n";
        }
        out << "\n";
    }
    else {
        std::cout << "Usage: estimateVarComp [SNP] [Pheno] [kinship] [outputPath]" << std::endl;
    }
    kinin.close();
    phenoin.close();
    snpin.close();
    out.close();
    return 0;
}

// read file to Eigen matrix
template<typename M>
M load_file(const std::string& path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    int rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ' ')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size() / rows);
}

// calc X'X
MatrixXd calc_xpx(const MatrixXd& X)
{
    const int n = X.cols();

    return MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(X.transpose());
}

// eigen decomposition
// returns eigenvalues and *transposed* eigenvectors
std::pair<VectorXd, MatrixXd> eigen_decomp(const MatrixXd& A)
{
    const SelfAdjointEigenSolver<MatrixXd> VLV(A);
    return std::make_pair(VLV.eigenvalues(), VLV.eigenvectors().transpose());
}

// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(const MatrixXd& K, const MatrixXd& y, const MatrixXd& X) {
    const std::pair<VectorXd, MatrixXd> e = eigen_decomp(K);
    const MatrixXd yrot = e.second * y;
    const MatrixXd Xrot = e.second * X;

    struct eigenrot result;
    result.Kva = e.first;
    result.Kve = e.second;
    result.y = yrot;
    result.X = Xrot;
    
    const std::pair<VectorXd, MatrixXd> ee = eigen_decomp(K);
    result.Kva = ee.first;
    MatrixXd Kve_t = ee.second;
    const MatrixXd yy = Kve_t * y;
    const MatrixXd XX = Kve_t * X;
    result.Kve = Kve_t;
    result.y = yy;
    result.X = XX;

    return result;
}

// calculate log det X'X
double calc_logdetXpX(const MatrixXd& X)
{
    const MatrixXd XpX(calc_xpx(X)); // calc X'X
    const int p = X.cols();

    // eigen decomposition of X'X
    const std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);

    // calculate log det X'X
    double result = 0.0;
    for (int i = 0; i < p; i++) result += log(e.first[i]);

    return result;
}

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
    const MatrixXd& X, const bool reml = true)
{
    const int n = Kva.size();
    const int p = X.cols();
    struct lmm_fit result;

    // diagonal matrix of weights
    VectorXd S(n);
    for (int i = 0; i < n; i++)
        S[i] = 1.0 / (hsq * Kva[i] + 1.0 - hsq);

    // calculate a bunch of matrices
    const MatrixXd XSt = X.transpose() * S.asDiagonal();
    MatrixXd ySt(1, n);
    for (int i = 0; i < n; i++) ySt(0, i) = y[i] * S[i];
    const MatrixXd XSX = XSt * X;
    const MatrixXd XSy = XSt * y;
    const MatrixXd ySy = ySt * y;

    // estimate of beta, by weighted LS
    const std::pair<VectorXd, MatrixXd>e = eigen_decomp(XSX);
    double logdetXSX = 0.0;
    VectorXd inv_evals(p);
    for (int i = 0; i < p; i++) {
        inv_evals[i] = 1.0 / e.first[i];
        if (reml) logdetXSX += log(e.first[i]);
    }
    const MatrixXd beta = e.second.transpose() * inv_evals.asDiagonal() * e.second * XSy;

    // residual sum of squares
    const MatrixXd rss = ySy - XSy.transpose() * beta;

    // return value
    result.rss = rss(0, 0);
    result.sigmasq = result.rss / (double)(n - p);
    result.beta = beta.col(0);
    result.logdetXSX = logdetXSX; // determinant (if REML)

    return result;
}

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
    const MatrixXd& X, const bool reml = true, const double logdetXpX = NULL)
{
    const int n = Kva.size();
    const int p = X.cols();

    // estimate beta and sigma^2
    struct lmm_fit ml_soln = getMLsoln(hsq, Kva, y, X, reml);

    // calculate log likelihood
    double loglik = (double)n * log(ml_soln.rss);
    for (int i = 0; i < n; i++)
        loglik += log(hsq * Kva[i] + 1.0 - hsq);
    loglik *= -0.5;

    if (reml) {
        double logdetXpX_val;
        if (logdetXpX == NULL) { // need to calculate it
            MatrixXd XpX(calc_xpx(X));
            std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);
            logdetXpX_val = 0.0;
            for (int i = 0; i < p; i++) logdetXpX_val += log(e.first[i]);
        }
        else logdetXpX_val = logdetXpX;

        loglik += 0.5 * (p * log(2 * M_PI * ml_soln.sigmasq) + logdetXpX_val - ml_soln.logdetXSX);
    }

    ml_soln.loglik = loglik;
    return ml_soln;
}

// just the negative log likelihood, for the optimization
double negLL(const double x, struct calcLL_args* args)
{
    const struct lmm_fit result = calcLL(x, args->Kva, args->y, args->X,
        args->reml, args->logdetXpX);

    return -result.loglik;
}

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
struct lmm_fit fitLMM(const VectorXd& Kva, const VectorXd& y, const MatrixXd& X, const bool reml = true, const bool check_boundary = true, 
                      const double logdetXpX = NULL, const double tol = 1e-4) {
    struct lmm_fit result;

    // calculate log det XpX, if necessary
    // (note same befor and after it's "rotated" by eigenvec of kinship matrix
    double logdetXpX_val;
    if (reml && logdetXpX == NULL) {
        MatrixXd XpX(calc_xpx(X));
        std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);
        int p = X.cols();
        logdetXpX_val = 0.0;
        for (int i = 0; i < p; i++) logdetXpX_val += log(e.first[i]);
    }
    else logdetXpX_val = logdetXpX;

    // function arguments for calcLL
    struct calcLL_args args;
    args.Kva = Kva;
    args.y = y;
    args.X = X;
    args.reml = reml;
    args.logdetXpX = logdetXpX_val;

    const double hsq = qtl2_Brent_fmin(0.0, 1.0, (double (*)(double, void*)) negLL, &args, tol);
    result = calcLL(hsq, Kva, y, X, reml, logdetXpX_val);
    result.hsq = hsq;

    if (check_boundary) {
        struct lmm_fit boundary_result;
        boundary_result = calcLL(0.0, Kva, y, X, reml, logdetXpX_val);
        if (boundary_result.loglik > result.loglik) {
            result = boundary_result;
            result.hsq = 0.0;
        }
        boundary_result = calcLL(1.0, Kva, y, X, reml, logdetXpX_val);
        if (boundary_result.loglik > result.loglik) {
            result = boundary_result;
            result.hsq = 1.0;
        }
    }

    return result;
}

double qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void*),
    void* info, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

    /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

    /*  main loop starts here ----------------------------------- */

    for (;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;

        /* check stopping criterion */

        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */

            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }

        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */

            d = p / q;
            u = x + d;

            /* f must not be evaluated too close to ax or bx */

            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }

        /* f must not be evaluated too close to x */

        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;

        fu = (*f)(u, info);

        /*  update  a, b, v, w, and x */

        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        }
        else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */

    return x;
}



int count_matrix_col(std::ifstream& matrix) {
    std::string token;
    std::stringstream stream; int count = 0;
    std::string read_line;

    std::getline(matrix, read_line);
    stream.str(read_line);
    while (stream >> token) {
        count++;
    }
    matrix.clear(); //coursor reset
    matrix.seekg(0, std::ios_base::beg);
    return count;
}

int count_matrix_row(std::ifstream& matrix) {
    std::string read_line; int count = 0;

    while (matrix.peek() != EOF) {
        std::getline(matrix, read_line);
        count++;
    }
    matrix.clear();//coursor reset
    matrix.seekg(0, std::ios_base::beg);
    return count;
}