
#include "MultiTrans.h"

//#include <algorithm>



int main() {

    std::cout << "@----------------------------------------------------------@" << std::endl;
    std::cout << "|       MTVCP      |      v1.0      |      10/12/2019      |"<< std::endl;
    std::cout << "| (C) 2019 Gi Ju Lee, Tae Gun Kim and Jong Wha Joanne Joo  |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| For documentation, citation & bug-report instructions:   |" << std::endl;
    std::cout << "|                                                          |" << std::endl;
    std::cout << "|             https://github.com/2eding/MTVCP              |" << std::endl;
    std::cout << "@----------------------------------------------------------@" << std::endl;

    //add error handling
    std::ifstream input_x("X.txt");
    std::ifstream input_y("Y.txt");
    std::ofstream output_c("result.txt");
    std::vector<double> vc_1;
    std::vector<double> vc_2;

    Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
    Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));
    Eigen::MatrixXd K = estimateKinship(X);
    estimateVarComp(K, X, Y, vc_1, vc_2); //calculate vc_1, vc_2
    Eigen::MatrixXd cor_matrix = cor_mat(K, X, vc_1, vc_2);
    //std::cout << VC <<std::endl;
    output_c << cor_matrix << std::endl;

    vc_1.~vector();
    vc_2.~vector();
    input_x.close();
    input_y.close();
    output_c.close();

    return 0;
}

Eigen::MatrixXd read_mat(std::ifstream& input_file, int row, int col) {

    std::string read_buffer;
    std::string token;
    std::stringstream stream;
    Eigen::MatrixXd ret_mat = Eigen::MatrixXd(row, col);

    for (int i = 0; i < row; i++) { //row
        std::getline(input_file, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < col; j++) { //col
            stream >> token;
            ret_mat(i, j) = std::stold(token);
        }
        stream.clear();  //for coursor reset
    }
    return ret_mat;
}

Eigen::MatrixXd estimateKinship(Eigen::MatrixXd snp) {
    std::ofstream out("Ktt.txt"); // output kinship file
    snp = snp.transpose().eval(); //    W.transposeInPlace(); (same)
    
    int n = snp.rows();
    int m = snp.cols();

    // Calculate mean
    Eigen::RowVectorXd mean = snp.colwise().mean();
    Eigen::MatrixXd matrixMean = mean;

    // Calculate variance
    Eigen::RowVectorXd var = (snp.rowwise() - mean).array().square().colwise().mean(); //
    Eigen::MatrixXd matrixVar = var;

    // Calculate standardDeviation
    Eigen::MatrixXd matrixStd = matrixVar.array().sqrt();
    Eigen::MatrixXd kinship;

    // Standardization
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            snp(j, i) = (snp(j, i) - matrixMean.col(i).value()) / matrixStd.col(i).value();
        }
    }
    // Estimate kinship
    // X * X^T / n
    kinship = snp * snp.transpose().eval() * 1.0 / m;
    out.precision(18);
    out << kinship << "\n";
    out.close();
    return kinship;
}
void estimateVarComp(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, Eigen::MatrixXd phenotype, std::vector<double>& vc_1, std::vector<double>& vc_2) {
    struct eigenrot e;
    struct lmm_fit vc;
    Eigen::MatrixXd ret_mat = Eigen::MatrixXd(phenotype.rows(), 2);
    for (int i = 0; i < phenotype.rows(); i++) {
        e = eigen_rotation(kinship, phenotype.row(i), snp);
        vc = fitLMM(e.Kva, e.y, e.X, true, true, NULL, 1e-6);
        vc_1.push_back(vc.hsq * vc.sigmasq);
        vc_2.push_back((1 - vc.hsq) * vc.sigmasq);
    }
}
Eigen::MatrixXd cor_mat(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, std::vector<double> vc_1, std::vector<double> vc_2) {
    double Median_Vg = cal_median(vc_1);
    double Median_Ve = cal_median(vc_2);
    Eigen::MatrixXd unit_mat;

    Eigen::MatrixXd ret(snp.rows(), snp.rows());
    std::vector<double> stand_devi;

    snp = snp.transpose().eval();
    kinship = (kinship * Median_Vg) + unit_mat.setIdentity(kinship.cols(), kinship.rows())*Median_Ve;
 
    kinship = kinship.sqrt();
    kinship = kinship.transpose().eval() * snp;
    cal_standard_deviation(stand_devi, kinship);
    cal_cor(stand_devi, kinship, ret);
    stand_devi.~vector();
    unit_mat.resize(0, 0);
    return ret;
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
double cal_median(std::vector<double> col) {
    std::vector<double> temp;
    double median = 0.0;
    for (int i = 0; i < col.size(); i++) {
        temp.push_back(col[i]);
    }
    std::sort(temp.begin(), temp.end());

    if (temp.size() % 2 == 0) { //even
        median = (temp[(temp.size() / 2)] + temp[(temp.size() / 2) - 1]) / 2;
    }
    else {
        median = temp[temp.size() / 2];
    }
    temp.~vector();
    return median;
}
void cal_standard_deviation(std::vector<double>& ret, Eigen::MatrixXd& mat) {
    double temp, temp_sum;
    for (int i = 0; i < mat.cols(); i++) { //col
        temp = 0;
        temp_sum = 0;
        for (int j = 0; j < mat.rows(); j++) { //row
            temp = mat(j, i) - (mat.col(i).sum() / ((double)mat.rows() - 1.0));
            temp_sum += temp*temp; // mat.col(i).sum()/(double)mat.rows() is col average
        }
        ret.push_back(sqrt((temp_sum)));
    }
}
void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat, Eigen::MatrixXd& colmat) {
    double temp, total_temp;
    for (int k = 0; k < mat.cols(); k++) { // cols
        for (int i = 0; i < mat.cols(); i++) { // cols
            total_temp = 0;
            for (int j = 0; j < mat.rows(); j++) { // rows 
                temp = mat(j, k) - (mat.col(k).sum() / (double)(mat.rows() - 1.0));
                temp = temp * (mat(j, i) - (mat.col(i).sum() / (double)(mat.rows() - 1.0)));
                total_temp += temp;
            }
            colmat(k, i) = (((total_temp))/(ret[k]*ret[i]));
        }
    }
}

// calc X'X
MatrixXd calc_xpx(const MatrixXd& X)
{
    const int n = X.cols();
    return MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(X.transpose().eval());
}

// eigen decomposition
// returns eigenvalues and *transposed* eigenvectors
std::pair<VectorXd, MatrixXd> eigen_decomp(const MatrixXd& A)
{
    const SelfAdjointEigenSolver<MatrixXd> VLV(A);
    return std::make_pair(VLV.eigenvalues(), VLV.eigenvectors().transpose().eval());
}

// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(const MatrixXd& K, const MatrixXd& y, const MatrixXd& X) {
    const std::pair<VectorXd, MatrixXd> e = eigen_decomp(K);
    // e.first's matrix size: # of individual x 1
    // e.second's matrix size: # of individual x # of individual
    const MatrixXd yrot = e.second * y.transpose().eval(); // y: 1 x individual
    const MatrixXd XX = MatrixXd(K.rows(), 1);
    const MatrixXd Xrot = e.second * XX; // X: snp x individual

    struct eigenrot result;
    result.Kva = e.first;
    result.Kve = e.second;
    result.y = yrot;
    result.X = Xrot;

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
    const MatrixXd XSt = X.transpose().eval() * S.asDiagonal();
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
    const MatrixXd beta = e.second.transpose().eval() * inv_evals.asDiagonal() * e.second * XSy;

    // residual sum of squares
    const MatrixXd rss = ySy - XSy.transpose().eval() * beta;

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
struct lmm_fit fitLMM(const VectorXd& Kva, const VectorXd& y, const MatrixXd& X, const bool reml = true, const bool check_boundary = true, const double logdetXpX = NULL, const double tol = 1e-4) {
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
double qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void*), void* info, double tol)
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

