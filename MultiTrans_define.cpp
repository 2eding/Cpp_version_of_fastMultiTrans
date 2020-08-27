#include "MultiTrans.h"


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

Eigen::MatrixXd estimateKinship(Eigen::MatrixXd X) {
    int snp_mat_row = X.rows();
    int snp_mat_col = X.cols();
    int n = snp_mat_col;
    int m = 1000;

    Eigen::MatrixXd W = Eigen::MatrixXd(n, m).setOnes() * NAN;
    int i = 0;
    Eigen::MatrixXd kinship = Eigen::MatrixXd(n, n).setZero();
    Eigen::MatrixXd kinship_j = Eigen::MatrixXd(n, n);
    while (i < snp_mat_row) {
        int j = 0;
        while ((j < m) && (i < snp_mat_row)) {
            double mean = X.row(i).mean();
            double snp_var = X.row(i).squaredNorm() - (mean * mean);
            if (snp_var == 0) {
                i += 1;
                continue;
            }
            W.col(j) = X.row(i).transpose();
            i += 1;
            j += 1;
        }
        if (j < m) {
            Eigen::MatrixXd W2 = Eigen::MatrixXd(n, j);
            for (int x = 0; x < n; x++) {
                for (int y = 0; y < j; y++) {
                    W2(x, y) = W(x, y);
                }
            }
            W.resize(n, j);
            W = W2;
        }
        if (kinship.isZero()) {
            kinship = calculateKinship(W.transpose()) * j;
        }
        else {
            kinship_j = calculateKinship(W.transpose()) * j;
            kinship = kinship + kinship_j;
        }
    }
    kinship = kinship / float(snp_mat_row);
    return kinship;
}

Eigen::MatrixXd calculateKinship(Eigen::MatrixXd W) {
    W = W.transpose().eval();
    int n = W.rows();
    int m = W.cols();

    // Calculate mean
    Eigen::RowVectorXd mean = W.colwise().mean();
    Eigen::MatrixXd matrixMean = mean;

    // Calculate variance
    Eigen::RowVectorXd var = (W.rowwise() - mean).array().square().colwise().mean();
    Eigen::MatrixXd matrixVar = var;

    // Calculate standardDeviation
    Eigen::MatrixXd matrixStd = matrixVar.array().sqrt();

    Eigen::MatrixXd kinship = Eigen::MatrixXd(m, m);

    // Standardization
    omp_set_num_threads(th_num);
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            W(j, i) = (W(j, i) - matrixMean.col(i).value()) / matrixStd.col(i).value();
        }
    }

    // Estimate kinship
    // X * X^T / n
    kinship = W * W.transpose().eval() * 1.0 / float(m);

    return kinship;
}

void estimateVarComp(Eigen::MatrixXd& kinship, Eigen::MatrixXd& snp, Eigen::MatrixXd& phenotype, std::vector<double>& vc_1, std::vector<double>& vc_2) {
    struct eigenrot e;
    struct lmm_fit vc;
    int pheno_rows = phenotype.rows();
    const std::pair<VectorXd, MatrixXd> ed = eigen_decomp(kinship);
    //omp_set_num_threads(th_num);
    //#pragma omp parallel for
    for (int i = 0; i < pheno_rows; i++) {
        e = eigen_rotation(kinship, phenotype.row(i), snp, ed);
        vc = fitLMM(e.Kva, e.y, e.X, true, true, 1e-6);
        vc_1.push_back(vc.hsq * vc.sigmasq);
        vc_2.push_back((1 - vc.hsq) * vc.sigmasq);
    }
    //    std::cout << cal_median(vc_1) << " " << cal_median(vc_2) << std::endl; // Vg Ve
}

void rotate_X_SigmaM(Eigen::MatrixXd& kinship, Eigen::MatrixXd& snp, double Vg, double Ve) {
    double Median_Vg = Vg;
    double Median_Ve = Ve;

    Eigen::MatrixXd unit_mat;
    std::vector<double> stand_devi;
    snp = snp.transpose().eval();
    kinship = (kinship * Median_Vg) + unit_mat.setIdentity(kinship.cols(), kinship.rows()) * Median_Ve;

    SelfAdjointEigenSolver<MatrixXd> es(kinship);

    kinship = es.operatorSqrt().inverse();
    kinship = kinship.transpose().eval() * snp;
}

//Eigen::MatrixXd cor_mat(Eigen::MatrixXd kinship, Eigen::MatrixXd snp, double Vg, double Ve) {
//
//    //double Median_Vg = Vg;
//    //double Median_Ve = Ve;
//    //
//    //Eigen::MatrixXd unit_mat;
//    //std::vector<double> stand_devi;
//    //std::cout << "snp transpose" << std::endl;
//    //snp = snp.transpose().eval();    
//    //std::cout << "kinship * Median_Vg and kinship.cols(), kinship.rows()" << std::endl;
//    //kinship = (kinship * Median_Vg) + unit_mat.setIdentity(kinship.cols(), kinship.rows()) * Median_Ve;
//    //
//    //SelfAdjointEigenSolver<MatrixXd> es(kinship);
//
//    //std::cout << "inverse matrix" << std::endl;
//
//    //kinship = es.operatorSqrt().inverse();
//    //kinship = kinship.transpose().eval() * snp;
//    //std::cout << "cal std" << std::endl;
//    //std::vector<double>* temp_std_vector;
//    //std::thread* th;
//    //for (int i = 0; i < 10; i++) {
//    //    std::cout << "thread i = " << i << std::endl;
//    //    if (i != 9) {
//    //        th[i] = std::thread(thread_func_std, i* (kinship.cols() / 10), (i + 1)* (kinship.cols() / 10), temp_std_vector[i], kinship);
//    //            //void thread_func(int begin, int end, std::vector<double> ret, Eigen::MatrixXd mat)
//    //    }
//    //    else {
//    //        th[i] = std::thread(thread_func_std, i* (kinship.cols() / 10), kinship.cols(), temp_std_vector[i], kinship);
//    //    }
//    //}
//    //std::cout << "thread start" << std::endl;
//    //for (int i = 0; i < 10; i++) {
//    //    th[i].join();
//    //}
//    //std::cout << "std_devi append" << std::endl;
//
//    //for (int i = 0; i < 10; i++) {
//    //    stand_devi.insert(stand_devi.end(), temp_std_vector[i].begin(), temp_std_vector[i].end());
//    //}
//    //std::cout << "std_devi write" << std::endl;
//    //std::ofstream st_dv("thread_stand_devi.txt");
//    //for (int i = 0; i < stand_devi.size(); i++) {
//    //    st_dv << stand_devi[i] << std::endl;
//    //}
//
//    //std::ifstream st_dv("stand_devi.txt");
//    //while (!st_dv.eof()) {
//    //    std::string token;
//    //    std::stringstream stream;
//    //    std::string r_b = "";
//    //    std::getline(st_dv, r_b);
//    //    if (r_b != "") {
//    //        stand_devi.push_back(std::stod(r_b)); //thread ¸ð»öÇØ¾ß´ï.
//    //    }
//    //}
//    //st_dv.close();
//
//    //unit_mat.resize(0, 0);
//    //snp.resize(0, 0);
//    //return cal_cor(stand_devi, kinship);
//
//}
void corrband_mat(int windowSize, std::ifstream& sig_cor, std::ofstream& corr_band) {
    std::string temp_line = "";
    std::string temp= "";
    std::stringstream ss;
    std::getline(sig_cor, temp_line); // first line delete;
    int snpCount = 0;
    while (!sig_cor.eof()) {
        snpCount++;
        std::getline(sig_cor, temp_line);
        ss.str(temp_line);
        if (temp_line != "") {
            for (int i = 0; i < (windowSize - snpCount); i++) {
                corr_band << "0 ";
            }
            for (int i = 0; i < snpCount - windowSize; i++) {
                ss >> temp;
            }
            if (snpCount <= windowSize) {
                for (int i = 0; i < snpCount; i++) { ss >> temp; corr_band << temp << " ";}
            }
            else {
                for (int i = 0; i < windowSize; i++) { ss >> temp; corr_band << temp << " "; }
            }
            if(!sig_cor.eof()){
            corr_band << "\n";
            }
        }
        ss.clear();
    }
    std::cout << snpCount << " number of lines read" << std::endl;
}
void thread_func_std(int begin, int end, int mat_rows, std::vector<double>& ret, Eigen::MatrixXd& mat){
    double temp, temp_sum, mat_col_sum_divide_kin_rows;
    int kin_rows = mat.rows();
    for (int i = begin; i < end; i++) { //col
        temp = 0;
        temp_sum = 0;
        mat_col_sum_divide_kin_rows = mat.col(i).sum() / ((double)kin_rows);
        //for (int j = 0; j < kin_rows; j++) { //row
        //    temp = mat(j, i) - (mat_col_sum_divide_kin_rows);
        //    temp_sum += temp * temp; // mat.col(i).sum()/(double)mat.rows() is col average
        //}
       temp_sum = mat.col(i).squaredNorm() - (mat_col_sum_divide_kin_rows * mat.col(i).sum());
       ret.push_back(sqrt((temp_sum)));
    }
    //std::cout << "vector size is " << ret.size() << std::endl;
}
void thread_func_std2(std::vector<double>& ret, std::vector<double>& mat_col_sum_divide_mat_rows, 
                                  int mat_cols, int mat_rows, int begin, int end, Eigen::MatrixXd& mat, int thread_number) {
    double total_temp;
    Eigen::MatrixXd temp;
    std::string path = "./multitrans_temporary_corrleationMatrix/cormat" + std::to_string(thread_number) + ".txt";
    std::ofstream divide_cor_mat(path);
    for (int k = begin; k < end; k++) { // cols
        temp = mat.col(k).transpose();
        for (int i = 0; i < mat_cols; i++) { // cols
            total_temp = 0;
            //for (int j = 0; j < mat_rows; j++) { // rows 
            //    temp = mat(j, k) - mat_col_sum_divide_mat_rows[k];
            //    temp = temp * (mat(j, i) - mat_col_sum_divide_mat_rows[i]);
            //    total_temp += temp;
            //    total_temp += mat(j, k) * mat(j, i);
            //}
            total_temp = (temp * mat.col(i)).sum();
            total_temp = total_temp - (mat_col_sum_divide_mat_rows[k] * (mat.col(i).sum())) - (mat_col_sum_divide_mat_rows[i] * (mat.col(k).sum())) + (mat_rows * mat_col_sum_divide_mat_rows[k] * mat_col_sum_divide_mat_rows[i]);
            //std::cout << "test k i " << k << " " << i << std::endl;
            divide_cor_mat.precision(16);
            if (i != mat_cols - 1) {
                divide_cor_mat << (((total_temp)) / (ret[k] * ret[i])) << " ";
            }
            else {
                divide_cor_mat << (((total_temp)) / (ret[k] * ret[i])) << std::endl;
            }
        }
    }
    divide_cor_mat.close();
}

std::vector<double> multi_cal_std2(MatrixXd& sigmaM,int sig_rows, int sig_cols) {
    std::vector<double> *temp_std_vector = new std::vector<double>[th_num];
    std::thread th[th_num];
    for (int i = 0; i < th_num; i++) {
        if (i != th_num-1) {
            th[i] = std::thread(thread_func_std, i* (sig_cols / th_num), (i + 1)* (sig_cols / th_num), sig_rows,std::ref(temp_std_vector[i]), std::ref(sigmaM));
                //void thread_func(int begin, int end, std::vector<double> ret, Eigen::MatrixXd mat)
        }
        else {
            th[i] = std::thread(thread_func_std, i* (sig_cols / th_num), sig_cols, sig_rows, std::ref(temp_std_vector[i]), std::ref(sigmaM));
        }
    }
    for (int i = 0; i < th_num; i++) {
        th[i].join();
    }
    for (int i = 1; i < th_num; i++) {
        temp_std_vector[0].insert(temp_std_vector[0].end(), temp_std_vector[i].begin(), temp_std_vector[i].end());
    }
    return temp_std_vector[0];
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

    if (temp.size() % 2 == 0 && temp.size() < 2) { //even
        median = (temp[(temp.size() / 2)] + temp[(temp.size() / 2) - 1]) / 2;
    }
    else {
        median = temp[temp.size() / 2];
    }
    return median;
}
//
//void cal_standard_deviation(Eigen::MatrixXd Sigma, Eigen::MatrixXd snp, double Vg, double Ve, std::ofstream& std_devi) {
//    double temp, temp_sum;
//    int Sigma_cols = Sigma.cols();
//    int Sigma_rows = Sigma.rows();
//    for (int i = 0; i < Sigma_cols; i++) { //col
//        temp = 0;
//        temp_sum = 0;
//        for (int j = 0; j < Sigma_rows; j++) { //row
//            temp = Sigma(j, i) - (Sigma.col(i).sum() / ((double)Sigma_rows));
//            temp_sum += temp * temp; // mat.col(i).sum()/(double)mat.rows() is col average
//        }
//        std_devi << sqrt(temp_sum) << std::endl;
//    }
//}
void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat) {
    double temp, total_temp;
    int mat_cols = mat.cols();
    int mat_rows = mat.rows();
    std::vector<double> mat_col_sum_divide_mat_rows;
    for (int i = 0; i < mat_cols; i++) {
        mat_col_sum_divide_mat_rows.push_back(mat.col(i).sum() / (double)(mat_rows));
    }
    // std::cout << "cal correlaton" << std::endl;
    //for (int k = 0; k < mat.cols(); k++) { // cols
    //    for (int i = 0; i < mat.cols(); i++) { // cols
    //        total_temp = 0;
    //        for (int j = 0; j < mat.rows(); j++) { // rows 
    //            temp = mat(j, k) - mat_col_sum_divide_mat_rows[k];
    //            temp = temp * (mat(j, i) - mat_col_sum_divide_mat_rows[i]);
    //            total_temp += temp;
    //        }
    //        colmat(k, i) = (((total_temp)) / (ret[k] * ret[i]));
    //        if (i % 1000 == 0) {
    //            std::cout << "i = " << i << std::endl;
    //        }
    //    }
    //    std::cout << " k = " << k << std::endl;
    //}
    std::thread th[th_num];
    for (int i = 0; i < th_num; i++) {
        if (i != th_num-1) {
            th[i] = std::thread(thread_func_std2, std::ref(ret), std::ref(mat_col_sum_divide_mat_rows), mat_cols, mat_rows, i * (mat_cols / th_num), (i + 1) * (mat_cols / th_num), std::ref(mat), i);
            //void thread_func(int begin, int end, std::vector<double> ret, Eigen::MatrixXd mat)
        }
        else {
            th[i] = std::thread(thread_func_std2, std::ref(ret), std::ref(mat_col_sum_divide_mat_rows), mat_cols, mat_rows, i * (mat_cols / th_num), mat_cols, std::ref(mat), i);
            //thread_func_std2(std::ref(ret), std::ref(mat_col_sum_divide_mat_rows), mat_cols, mat_rows, i * (mat_cols / 10), mat_cols, std::ref(mat), std::ref(colmat));
        }
    }
    std::cout << "thread start" << std::endl;
    for (int i = 0; i < th_num; i++) {
        th[i].join();
    }
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
struct eigenrot eigen_rotation(const MatrixXd& K, const MatrixXd& y, const MatrixXd& X, const std::pair<VectorXd, MatrixXd>& e) {
    // e.first's matrix size: # of individual x 1
    // e.second's matrix size: # of individual x # of individual
    const MatrixXd yrot = e.second * y.transpose(); // y: 1 x individual
    const MatrixXd XX = MatrixXd::Ones(K.rows(), 1);
    const MatrixXd Xrot = e.second * XX; // X: snp x individual

    struct eigenrot result;
    result.Kva = e.first;
    result.Kve = e.second;
    result.y = yrot;
    result.X = Xrot;

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
    result.sigmasq = result.rss / (double(n) - double(p));
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
    const MatrixXd& X, const bool reml = true, const double logdetXpX = 0.0)
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
        if (logdetXpX == 0.0) { // need to calculate it
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
struct lmm_fit fitLMM(const VectorXd& Kva, const VectorXd& y, const MatrixXd& X, const bool reml = true, const bool check_boundary = true, const double tol = 1e-4) {
    struct lmm_fit result;

    // calculate log det XpX, if necessary
    // (note same befor and after it's "rotated" by eigenvec of kinship matrix
    double logdetXpX_val;
    MatrixXd XpX(calc_xpx(X));
    std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);
    int p = X.cols();
    logdetXpX_val = 0.0;
    for (int i = 0; i < p; i++) logdetXpX_val += log(e.first[i]);

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
double qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void*), void* info, double tol){
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

FILE* readfile(char* filename) {
    FILE* fp;
    if ((fp = fopen(filename, "r")) == NULL) { fprintf(stderr, "Cannot read file %s\n", filename); exit(-1); }
    return(fp);
}
FILE* writefile(char* filename) {
    FILE* fp;
    if ((fp = fopen(filename, "w")) == NULL) { fprintf(stderr, "Cannot write to file %s\n", filename); exit(-1); }
    return(fp);
}


