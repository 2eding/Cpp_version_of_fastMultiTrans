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
Eigen::MatrixXd corrband_mat(int windowSize, Eigen::MatrixXd cor_mat) {
    Eigen::MatrixXd ret_mat(cor_mat.rows() - 1, windowSize);
    ret_mat.setZero();
    for (int i = 0; i < ret_mat.rows(); i++) {
        for (int j = 0; j < i+1; j++) {
            ret_mat(i, ret_mat.cols() - 1 + j - i) = cor_mat(i + 1, j);
        }
    }
    return ret_mat;
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


struct tables { double stat, prob; };

void slide_1prep(int windowsize, char* input_c, char* output)
{
    MKL_Set_Num_Threads(1); // avoid multi-threading
    int i, j, k, l, n11, n12, n21, n22, c_p, c_n, r0, r1, r2, g0, g1, g2;
    int nsnps, ncont, ncase, nindv, nmin, szwin, ws;
    int* a0, * a1, * a2, ntags, * tags, * connsec;
    char* bins, * ucpi, * ucpj;
    int lbound, rbound;
    char* binf, * outf, buf[1000], * line, * token;
    double f, * rs, * rsq, * convec, * constd, * scale, * thres, ** conthres, ** conscale, * denoms;
    double* vecx, * vecy, effntags;
    double cutpv, numer, max_st, * h1, * logfact, con, logfactcon, subt, stat, prob, statconst, cutstat, midp;
    double statcutconst_p, statcutconst_n, subcon, nextp, prevprob, maxprob;
    double negzcut, ltailp, * maxchi2, * negzs, * pointwiseps;
    double x = DOSE_EFF, xx = x * x; // dose effect
    int ione = 1, info, nrow, ncol;
    FILE* fp;
    struct tables* t_p, * t_n;
    enum MODE { H, G, QH, QG, C } mode;
    // receive arguments
    //if (argc < 2) print_usage();
    //if (!strcmp(argv[1], "-H")) { mode = H; }
    //else if (!strcmp(argv[1], "-G")) { mode = G; }
    //else if (!strcmp(argv[1], "-QH")) { mode = QH; }
    //else if (!strcmp(argv[1], "-QG")) { mode = QG; }
    //else if (!strcmp(argv[1], "-C")) { mode = C; }
    //else { print_usage(); }
    //if (mode == H || mode == G) {
    //  if (argc != 7) print_usage();
    //  binf = argv[2];
    //  ncont = atoi(argv[3]); // # controls
    //  ncase = atoi(argv[4]); // # cases
    //  nindv = ncont + ncase; 
    //  nmin = min(ncont, ncase);
    //  szwin = atoi(argv[5]); // window size
    //  outf = argv[6];
    //} else if (mode == QH || mode == QG) {
    //  if (argc != 6) print_usage();
    //  binf = argv[2];
    //  nindv = atoi(argv[3]);
    //  szwin = atoi(argv[4]); // window size
    //  outf = argv[5];
    //  ncase = ncont = nmin = 0;
    //} else if (mode == C) {
    //  if (argc != 5) print_usage();
    //  binf = argv[2];
    //  szwin = atoi(argv[3]);
    //  outf = argv[4];
    //  ncase = ncont = nmin = nindv = 0;
    //}
    //test 2020/04/01
    mode = C;
    binf = input_c; // input c.txt
    szwin = windowsize; //windowsize
    outf = output; //output file ex) prep
    ncase = ncont = nmin = nindv = 0;
    //
    //count # of row and column (and sanity check)
    fprintf(stderr,">>> Counting SNPs\n");
    nrow = 0;
    line = (char*)malloc(20000001 * sizeof(char));
    fp = readfile(binf);
    while (fgets(line, 20000000, fp) != NULL) {
        nrow++;
        token = strtok(line, " \t\n");
        if (token == NULL || token[0] == '#') nrow--;
    }
    if (nrow == 0) { fprintf(stderr, "File empty\n"); exit(-1); }
    fclose(fp);
    nsnps = nrow;
    if (mode == C) {
        nsnps++;
        ntags = nsnps;
    }
    fprintf(stderr, ">>> Number of SNPs: %d\n", nsnps);
    // memory alloc
    bins = (char*)malloc(sizeof(char) * (szwin + 1) * nindv);
    tags = (int*)malloc(sizeof(int) * nsnps);
    a0 = (int*)malloc(sizeof(int) * nsnps);
    a1 = (int*)malloc(sizeof(int) * nsnps);
    a2 = (int*)malloc(sizeof(int) * nsnps);
    h1 = (double*)malloc(sizeof(double) * nsnps);
    denoms = (double*)malloc(sizeof(double) * nsnps);
    rs = (double*)calloc(sizeof(double), szwin * szwin);
    rsq = (double*)malloc(sizeof(double) * (szwin) * (szwin));//denotes square matrix of r.. not r^2.
    vecx = (double*)malloc(sizeof(double) * (szwin));
    vecy = (double*)malloc(sizeof(double) * (szwin));
    constd = (double*)malloc(sizeof(double) * nsnps);
    convec = (double*)malloc(sizeof(double) * nsnps * szwin);
    conscale = (double**)malloc(sizeof(double*) * 2 * nsnps);
    conthres = (double**)malloc(sizeof(double*) * 2 * nsnps);
    scale = (double*)malloc(sizeof(double) * MAXNSEC);
    thres = (double*)malloc(sizeof(double) * MAXNSEC);
    connsec = (int*)malloc(sizeof(int) * 2 * nsnps);
    t_p = (struct tables*)malloc((nmin + 2) * (nmin + 1) / 2 * sizeof(struct tables));
    t_n = (struct tables*)malloc((nmin + 2) * (nmin + 1) / 2 * sizeof(struct tables));
    logfact = (double*)malloc(sizeof(double) * (nindv + 1));
    negzs = (double*)malloc(sizeof(double) * nsnps);
    pointwiseps = (double*)malloc(sizeof(double) * nsnps);
    maxchi2 = (double*)malloc(sizeof(double) * nsnps);

    if (mode == C) {
        for (i = 0; i < ntags; ++i) {  // unit scaling.
            l = 0;
            thres[l] = -100;
            scale[l] = 1;
            l++;
            connsec[2 * i] = l;
            conthres[2 * i] = (double*)malloc(sizeof(double) * (l));
            conscale[2 * i] = (double*)malloc(sizeof(double) * (l));
            memcpy(conthres[2 * i], thres, sizeof(double) * l);
            memcpy(conscale[2 * i], scale, sizeof(double) * l);
            l = 0;
            thres[l] = 100;
            scale[l] = 1;
            l++;
            connsec[2 * i + 1] = l;
            conthres[2 * i + 1] = (double*)malloc(sizeof(double) * (l));
            conscale[2 * i + 1] = (double*)malloc(sizeof(double) * (l));
            memcpy(conthres[2 * i + 1], thres, sizeof(double) * l);
            memcpy(conscale[2 * i + 1], scale, sizeof(double) * l);
        }
    }// end of scaling factor part

    //Computing, finally, convec and constd.
    fprintf(stderr, ">>> Computing conditional distributions\n");
    constd[0] = 1.0;
    fp = readfile(binf); // read binary file (2nd time)

    for (i = 1; i < ntags; ++i) {// 'i' is of interest
        if (i > 1 && i % 10000 == 0) fprintf(stderr, "SNP %d processed\n", i);
        ws = min(szwin, i);//effective window lookup size
        // compute "r" for one more line
        
        if (mode == C) {
            for (j = 0; j < szwin; ++j) {
                fscanf(fp, "%lf", &rs[(i % szwin) * szwin + j]);
            }
        }
        for (j = 0; j < ws; ++j) {
            rsq[j * ws + j] = 1.0 + EIGEN_EPS;
            for (k = 0; k < j; ++k) { rsq[j * ws + k] = rs[((i - ws + j) % szwin) * szwin + (szwin - j + k)]; }
            vecx[j] = vecy[j] = rs[(i % szwin) * szwin + szwin - ws + j];
        }
        dpotrf("U", &ws, rsq, &ws, &info);
        dpotrs("U", &ws, &ione, rsq, &ws, vecy, &ws, &info);
        f = ddot(&ws, vecx, &ione, vecy, &ione);
        memset(&convec[i * szwin], 0, sizeof(double) * (szwin - ws));
        memcpy(&convec[i * szwin + szwin - ws], vecy, sizeof(double) * ws);
        constd[i] = sqrt(max(1. + EIGEN_EPS - f, 0));
        //printf("%f\n", constd[i]);
    }
    fclose(fp);

    // Write the data.
    sprintf(buf, "%s.txt_polySNPs", outf);
    fp = writefile(buf);
    for (i = 0; i < ntags; ++i) fprintf(fp, "%d\n", tags[i]);
    fclose(fp);

    sprintf(buf, "%s.bin_constd", outf);
    fp = writefile(buf);
    fwrite(constd, sizeof(double), ntags, fp);
    fclose(fp);

    sprintf(buf, "%s.bin_convec", outf);
    fp = writefile(buf);
    fwrite(convec, sizeof(double), ntags * szwin, fp);
    fclose(fp);

    /// print purpose
 /*  sprintf(buf,"%s.constd.txt",outf); 
     fp = writefile(buf); 
     for(i=0;i < ntags;++i) fprintf(fp, "%.3lf\n", constd[i]); 
     fclose(fp); 
     sprintf(buf,"%s.convec.txt",outf); 
     fp = writefile(buf); 
     for(i=0;i < ntags;++i) { 
       for(k=0; k < szwin;k++) { fprintf(fp, "%.3lf ", convec[i*szwin+k]); } 
       fprintf(fp, "\n"); 
     } 
    fclose(fp); */

    sprintf(buf, "%s.txt_nsec", outf); //"%s.txt_nsec"
    fp = writefile(buf);
    for (i = 0; i < 2 * ntags; ++i) fprintf(fp, "%d\n", connsec[i]);
    fclose(fp);

    sprintf(buf, "%s.bin_scale", outf); //"%s.bin_scale"
    fp = writefile(buf);
    for (i = 0; i < 2 * ntags; ++i) {
        fwrite(&conthres[i][0], sizeof(double), connsec[i], fp);
        fwrite(&conscale[i][0], sizeof(double), connsec[i], fp);
    }
    fclose(fp);

    sprintf(buf, "%s.txt_wsize", outf); //"%s.txt_wsize"
    fp = writefile(buf);
    fprintf(fp, "%d\n", szwin);
    fclose(fp);

    if (mode == H || mode == G) { // each SNP's statistic
      //    vmlSetMode(VML_HA);
        vdCdfNorm(nsnps, negzs, pointwiseps);
        sprintf(buf, "%s.txt_pointwisep", outf);
        fp = writefile(buf);
        //    for(i=0;i < nsnps;++i) fprintf(fp, "%.10lf\t%.10lf\n", negzs[i]*negzs[i], pointwiseps[i]*2);
        for (i = 0; i < nsnps; ++i) {
            fprintf(fp, "%.10le\n", pointwiseps[i] * 2);
            //      fprintf(stdout, "%.10le %.10le\n", negzs[i], pointwiseps[i]*2);
        }
        fclose(fp);
    }

    sprintf(buf, "%s.txt_log", outf);
    fp = writefile(buf);
    fprintf(fp, ">>> Command invoked:");
    //for (i = 0;i < argc;i++) fprintf(fp," %s",argv[i]); fprintf(fp,"\n");
    fprintf(fp, ">>> Input file %s has been successfully pre-processed\n", binf);
    fprintf(fp, ">>> Number of SNPs: %d\n", nsnps);
    fprintf(fp, ">>> Number of removed non-polymorphic SNPs: %d\n", nsnps - ntags);
    fprintf(fp, ">>> Number of polymorphic SNPs processed: %d\n", ntags);
    fprintf(fp, ">>> The following files were created.\n");
    fprintf(fp, ">>> %s.bin_convec\t:\tconditional distribution information\n", outf);
    fprintf(fp, ">>> %s.bin_constd\t:\tconditional distribution information\n", outf);
    fprintf(fp, ">>> %s.bin_scale\t:\tscaling factor information\n", outf);
    fprintf(fp, ">>> %s.txt_nsec\t:\tscaling factor information\n", outf);
    fprintf(fp, ">>> %s.txt_polySNPs\t:\tindex of polymorphic SNPs\n", outf);
    fprintf(fp, ">>> %s.txt_wsize\t:\tstore window size\n", outf);
    if (mode == H || mode == G) fprintf(fp, ">>> %s.txt_pointwisep\t:\tpointwise p-values of every SNP\n", outf);
    fprintf(fp, ">>> %s.txt_log\t:\tlog file of this preprocess\n", outf);
    fclose(fp);

    fprintf(stdout, ">>> Input file %s has been successfully pre-processed\n", binf);
    fprintf(stdout, ">>> Number of removed non-polymorphic SNPs: %d\n", nsnps - ntags);
    fprintf(stdout, ">>> Number of polymorphic SNPs processed: %d\n", ntags);
    fprintf(stdout, ">>> The following files were created.\n");
    fprintf(stdout, ">>> %s.bin_convec\t:\tconditional distribution information\n", outf);
    fprintf(stdout, ">>> %s.bin_constd\t:\tconditional distribution information\n", outf);
    fprintf(stdout, ">>> %s.bin_scale\t:\tscaling factor information\n", outf);
    fprintf(stdout, ">>> %s.txt_nsec\t:\tscaling factor information\n", outf);
    fprintf(stdout, ">>> %s.txt_polySNPs\t:\tindex of polymorphic SNPs\n", outf);
    fprintf(stdout, ">>> %s.txt_wsize\t:\tstore window size\n", outf);
    if (mode == H || mode == G) fprintf(stdout, ">>> %s.txt_pointwisep\t:\tpointwise p-values of every SNP\n", outf);
    fprintf(stdout, ">>> %s.txt_log\t:\tlog file of this preprocess\n", outf);
}

int compare_by_stat(const void* a, const void* b) {
    double temp = ((struct tables*)a)->stat - ((struct tables*)b)->stat;
    return (temp > 0.) ? (1) : ((temp < 0.) ? (-1) : 0);
}

void slide_2run(char* input, char* output)
{
    MKL_Set_Num_Threads(1); // avoid multi-threading
    int ntags, szwin, nrand, * nsec;
    int i, j, l, seed, nloop, nsets, set, loop;
    char* inf, * outf;
    char buf[1000], * line;
    FILE* fp;
    double r, m, f;
    double* convec, * constd, ** scale, ** thres, * rnds, * maxnull, * rgen, * rconstd;
    VSLStreamStatePtr stream;
    int ione = 1;

    inf = input;
    outf = output;
    nrand = 100; //need input
    seed = 123;

    // count # of SNPs
    ntags = 0;
    line = (char*)malloc(20000001 * sizeof(char));
    sprintf(buf, "%s.txt_polySNPs", inf);
    fp = readfile(buf);
    while (fgets(line, 20000000, fp) != NULL) ntags++;
    fclose(fp);
    // read window size
    sprintf(buf, "%s.txt_wsize", inf);
    fp = readfile(buf);
    fscanf(fp, "%d", &szwin);
    fclose(fp);
    // alloc memory
    convec = (double*)malloc(sizeof(double) * ntags * szwin);
    constd = (double*)malloc(sizeof(double) * ntags);
    thres = (double**)malloc(sizeof(double*) * 2 * ntags);
    scale = (double**)malloc(sizeof(double*) * 2 * ntags);
    nsec = (int*)malloc(sizeof(int) * 2 * ntags);
    // read binaries
    sprintf(buf, "%s.bin_convec", inf);
    fp = readfile(buf);
    fread(convec, sizeof(double), ntags * szwin, fp);
    fclose(fp);
    sprintf(buf, "%s.bin_constd", inf);
    fp = readfile(buf);
    fread(constd, sizeof(double), ntags, fp);
    fclose(fp);
    sprintf(buf, "%s.txt_nsec", inf);
    fp = readfile(buf);
    for (i = 0; i < 2 * ntags; i++) fscanf(fp, "%d", &nsec[i]);
    fclose(fp);
    sprintf(buf, "%s.bin_scale", inf);
    fp = readfile(buf);
    for (i = 0; i < 2 * ntags; i++) {
        thres[i] = (double*)malloc(sizeof(double) * nsec[i]);
        scale[i] = (double*)malloc(sizeof(double) * nsec[i]);
        fread(&thres[i][0], sizeof(double), nsec[i], fp);
        fread(&scale[i][0], sizeof(double), nsec[i], fp);
    }
    fclose(fp);
    // prepare for sampling
    vslNewStream(&stream, VSL_BRNG_MT19937, seed);
    rnds = (double*)malloc(sizeof(double) * ntags);
    rgen = (double*)malloc(sizeof(double) * ntags);
    rconstd = (double*)malloc(sizeof(double) * ntags);
    maxnull = (double*)malloc(sizeof(double) * ONEWRITE * 3);

    fp = writefile(outf);
    nsets = max(1, nrand / ONEWRITE);
    for (set = 0; set < nsets; set++) {
        if (set == nsets - 1) {
            nloop = nrand - set * ONEWRITE;
        }
        else {
            nloop = ONEWRITE;
        }
        if (set > 0) fprintf(stderr, "Sampling %d\n", set * ONEWRITE);
        for (loop = 0; loop < nloop; ++loop) {
            m = 0;
            //FIX VSL_METHOD_DGAUSSIAN_BOXMULLER -> VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
            vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, ntags, rgen, 0., 1.);
            vdMul(ntags, rgen, constd, rconstd);
            for (j = 0; j < ntags; ++j) {
                r = rconstd[j];
                if (j < szwin) {
                    r += ddot(&j, &convec[j * szwin + szwin - j], &ione, &rnds[0], &ione);
                   // printf("test2, %f \n", r);
                }
                else {
                    r += ddot(&szwin, &convec[j * szwin], &ione, &rnds[j - szwin], &ione);
                    //printf("test2-2, %f \n", r);
                }
                rnds[j] = r;
                f = r;
                if (f <= 0) { // negative case
                    if (nsec[2 * j] > 0) {
                        if (f >= thres[2 * j][0]) { // if 1st section.
                            f *= scale[2 * j][0];
                        }
                        else {
                            for (l = 1; l < nsec[2 * j]; l++) {
                                if (f > thres[2 * j][l]) {
                                    f *= scale[2 * j][l - 1] + (scale[2 * j][l] - scale[2 * j][l - 1]) / (thres[2 * j][l] - thres[2 * j][l - 1]) * (f - thres[2 * j][l - 1]); // linear interpolation
                                    break;
                                }
                            }
                            if (l == nsec[2 * j]) f *= scale[2 * j][l - 1]; // over the final section
                        }
                    }
                }
                else { // positive case
                    if (nsec[2 * j + 1] > 0) {
                        if (f <= thres[2 * j + 1][0]) { // if 1st section
                            f *= scale[2 * j + 1][0];
                        }
                        else {
                            for (l = 1; l < nsec[2 * j + 1]; l++) {
                                if (f < thres[2 * j + 1][l]) {
                                    f *= scale[2 * j + 1][l - 1] + (scale[2 * j + 1][l] - scale[2 * j + 1][l - 1]) / (thres[2 * j + 1][l] - thres[2 * j + 1][l - 1]) * (f - thres[2 * j + 1][l - 1]); // linear interpolation
                                    break;
                                }
                            }
                            if (l == nsec[2 * j + 1]) f *= scale[2 * j + 1][l - 1]; // over the final section
                        }
                    }
                }
                if (fabs(f) > m && _finite(fabs(f))) { m = fabs(f); } //....&& _finite(fabs(f))
                else if ( m == 0 && !_finite(fabs(f))) { m = 0; }
            } //for(j)
            maxnull[loop] = m;
        } //for(loop)
        fwrite(maxnull, sizeof(double), nloop, fp);
             //for(i=0;i<nloop;i++)printf("test plz.. : %.10lf\n",maxnull[i]);
    }//for(set)
    fclose(fp);
    fprintf(stdout, ">>> Sampling result file %s has been successfully created.\n", outf);
    free(convec);
    free(constd);
    free(thres);
    free(scale);
    free(rnds);
    free(rgen);
    free(rconstd);
    free(maxnull);
}


void slide_3sort(char* input_max, char* output_sor)
{
    MKL_Set_Num_Threads(1); // avoid multi-threading
    char* orgf;
    int i, j, size, prevsize;
    FILE* fp;
    double chunk[DBL_CHUNK], * maxstats, * stats;


    orgf = output_sor; //output

    // # sampling consistency check
    prevsize = -1;
    fp = readfile(input_max);
    size = 0;
    do {
        j = fread(chunk, sizeof(double), DBL_CHUNK, fp);
        size += j;
    } while (j == DBL_CHUNK);
    if (prevsize != -1 && prevsize != size) {
        fprintf(stderr, "Number of sampling has to be the same for all result files\n"); exit(-1);
    }
    prevsize = size;
    fclose(fp);


    // read result files
    maxstats = (double*)malloc(sizeof(double) * size);
    stats = (double*)malloc(sizeof(double) * size);
    for (j = 0; j < size; j++) maxstats[j] = 0.;
    for (i = 2; i < 3; i++) {
        fp = readfile(input_max);
        fread(stats, sizeof(double), size, fp);
        fclose(fp);
        for (j = 0; j < size; j++) {
            maxstats[j] = max(maxstats[j], stats[j]);
            printf("test_max=%d: %f \n", j, maxstats[j]);
        }
    }

    // sort
    qsort(maxstats, size, sizeof(double), compare_double);

    // write 
    fp = writefile(orgf);
    fwrite(maxstats, sizeof(double), size, fp);
    fclose(fp);
}

int compare_double_rev(const void* s, const void* t) {
    double a = *((double*)s);
    double b = *((double*)t);
    return (a > b) ? -1 : a < b ? 1 : 0;
}
int compare_double(const void* s, const void* t) {
    double a = *((double*)s);
    double b = *((double*)t);
    return (a > b) ? 1 : a < b ? -1 : 0;
}

void slide_4correct(char* input_sor, char* output)
{
    MKL_Set_Num_Threads(1); // avoid multi-threading
    char* orgf, * outf, * mapf, * dataf;
    int i, j, size, nrow, ndata;
    FILE* fp;
    double chunk[DBL_CHUNK], * maxstats, x;
    double* data, * result1, * result2, dat, resul1, resul2, neff;
    char* name, * line, * token;

    enum MODE { P, T, IP, IT } mode;
    mode = P;
    // receive arguments
    //if (argc < 3) print_usage();
    //if (!strcmp(argv[1], "-p")) { mode = P; }
    //else if (!strcmp(argv[1], "-t")) { mode = T; }
    //else if (!strcmp(argv[1], "-ip")) { mode = IP; }
    //else if (!strcmp(argv[1], "-it")) { mode = IT; }
    //else { print_usage(); }
    orgf = input_sor;  //input sorted file
    // read orgf
    fp = readfile(orgf);
    size = 0;
    do {
        j = fread(chunk, sizeof(double), DBL_CHUNK, fp);
        size += j;
    } while (j == DBL_CHUNK);
    fclose(fp);
    if (size == 0) { fprintf(stderr, "orgainzed file empty\n"); exit(-1); }
    maxstats = (double*)malloc(sizeof(double) * size);
    fp = readfile(orgf);
    fread(maxstats, sizeof(double), size, fp);
    fclose(fp);

    // by modes
    dataf = "threshold.txt"; //threshold
    // read data
    fp = readfile(dataf);
    ndata = 0;
    while (fscanf(fp, "%lf", &x) != EOF) ndata++;
    fclose(fp);
    data = (double*)malloc(sizeof(double) * ndata);
    result1 = (double*)malloc(sizeof(double) * ndata);
    result2 = (double*)malloc(sizeof(double) * ndata);
    fp = readfile(dataf);
    for (i = 0; i < ndata; ++i) fscanf(fp, "%lf", &data[i]);
    fclose(fp);
    outf = output; //outputfile
    //for (i = 0;i < 10;i++) fprintf(stderr,"test_final? : %.10lf\n",data[i]);
    // map file 

    fp = writefile(outf);
    fprintf(fp, "#Result based on %d sampling\n", size);
    if (mode == P) {
        for (i = 0; i < ndata; i++) {
            correct_pvalue(maxstats, size, data[i], &result1[i], &result2[i]);
            /*printf("test re_1 : %f \n", result1[i]);
            printf("test re_2 : %f \n", result2[i]);*/
        }
        fprintf(fp, "#SNP_id\tPointwise-P\tCorrected-P\tApprox.STDEV(%% of P)\tEffective#ofTests\tNote\n");
    }
    for (i = 0; i < ndata; i++) {
        if (mode == P) {
            fprintf(fp, "SNP%d", i + 1);
        }
        fprintf(fp, "\t%.8le", data[i]);
        fprintf(fp, "\t%.8le", result1[i]);
        fprintf(fp, "\t%.8le", result2[i]);
        neff = (mode == P) ? result1[i] / data[i] : data[i] / result1[i];
        if (result1[i] < EPS || result1[i] >= 1. - EPS) {
            fprintf(fp, "(%.5lf%%)", 0.);
            fprintf(fp, "\t----");
        }
        else {
            fprintf(fp, "(%.5lf%%)", result2[i] / result1[i] * 100);
            fprintf(fp, "\t%.0lf", neff);
        }
        if (result1[i] < EPS) {
            fprintf(fp, "\tNeed more sampling");
        }
        else if (result1[i] >= 1. - EPS) {
            fprintf(fp, "\t");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

}

void correct_pvalue(double* maxstats, int size, double data, double* result1, double* result2) {
    int lfin, rfin, mid, x;
    double ltailp, negz, absz;
    ltailp = data / 2.;
    vdCdfNormInv(1, &ltailp, &negz);
    absz = (-1.) * negz;
    lfin = 0;
    rfin = size - 1;
    do {
        mid = (lfin + rfin) / 2;
        if (absz < maxstats[mid]) {
            rfin = mid;
        }
        else if (absz > maxstats[mid]) {
            lfin = mid;
        }
        else if (fabs(absz - maxstats[mid]) < EPS) {
            break;
        }
    } while (rfin - lfin > 2);
    x = mid;

    if (absz <= maxstats[x]) {
        while (x >= 0 && absz <= maxstats[x]) --x;
    }
    else if (absz > maxstats[x]) {
        while (x < size && absz > maxstats[x]) ++x;
        if (x < size) --x;
    }

    if (x == size) { *result1 = 0.; }
    else { *result1 = ((double)size - (x + 1)) / size; }
    *result2 = sqrt((*result1) * (1 - (*result1)) / size);
    return;
}

void per_marker_threshold(double* maxstats, int size, double data, double* result1, double* result2) {
    double negz, ltailp;
    negz = (-1.) * maxstats[size - (int)(size * data)];
    vdCdfNorm(1, &negz, &ltailp);
    *result1 = 2. * ltailp;
    *result2 = sqrt((*result1) * (1 - (*result1)) / size);
    return;
}
