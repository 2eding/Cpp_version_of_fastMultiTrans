#include "MultiTrans.h"

//#include <algorithm>

int main(int argc, char** argv) {

    std::cout << "@----------------------------------------------------------@" << std::endl;
    std::cout << "|       MTVCP      |      v1.0      |      10/12/2019      |" << std::endl;
    std::cout << "| (C) 2019 Gi Ju Lee, Tae Gun Kim and Jong Wha Joanne Joo  |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "|                                                          |" << std::endl;
    std::cout << "| For documentation, citation & bug-report instructions:   |" << std::endl;
    std::cout << "|             https://github.com/2eding/MTVCP              |" << std::endl;
    std::cout << "@----------------------------------------------------------@" << std::endl;

    //add error handling
   
    std::ifstream input_x(argv[1]);
    std::ifstream input_y(argv[2]);
    std::ofstream output_corrband(argv[3]);
    std::string r_path = argv[4];
    int windowsize = std::stoi(argv[5]);
    std::vector<double> vc_1;
    std::vector<double> vc_2;
    double median_vc1 = 0.0;
    double median_vc2 = 0.0;
    int sig_rows = 0;
    int sig_cols = 0;
    
    std::vector<double> std_vector;
    std::string command = "cat ";

    if (argc != 6) {
        std::cout << "please input genotypeData, PhenotypeData, corrbandmatrix_outputpath, correlation_output_path, windowsize, " << std::endl;
        return 0;
    }
    //std::cout << "input x y" << std::endl;

    Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
    Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));

    std::cout << "Estimate Kinship" << std::endl;
    Eigen::MatrixXd K = estimateKinship(X);
    std::cout << "EstimateVarComp" << std::endl;
    estimateVarComp(K, X, Y, vc_1, vc_2); //calculate vc_1, vc_2
    std::cout << "Calculate Median" << std::endl;
    median_vc1 = cal_median(vc_1); median_vc2 = cal_median(vc_2);

    //std::ofstream temp_vc("c++version_vc.txt");
    //temp_vc.precision(16);
    //for (int i = 0; i < vc_1.size(); i++) {
    //    temp_vc << vc_1[i] << " " << vc_2[i] << std::endl;
    //}
    //temp_vc.close();

    std::cout << "Calculate SigmaMatrix" << std::endl;
    rotate_X_SigmaM(K, X, median_vc1, median_vc2); // after this function, K is sigmaMatrix
    sig_rows = K.rows(); sig_cols = K.cols();
    
    std::cout << "Calculate standard deviation" << std::endl;
    std_vector = multi_cal_std2(K, sig_rows, sig_cols);
    std::cout << "Calculate Correlation Matrix" << std::endl;
    std::system("mkdir ./multitrans_temporary_corrleationMatrix");
    
    cal_cor(std_vector, K);

    for (int i = 0; i < th_num; i++) {
        command += "./multitrans_temporary_corrleationMatrix/cormat" + std::to_string(i) + ".txt ";
    }
    command += ">> " + r_path;
    std::system(command.c_str()); //merge correlation file command;
    std::system("rm -rf ./multitrans_temporary_corrleationMatrix");

    std::ifstream input_sigma_cor(r_path);
    if (!input_sigma_cor.is_open()) {
        std::cout << "cannot find correlation matrix file r.txt" << std::endl; return 0;
    }
    std::cout << "Calculate Corrband Matrix" << std::endl;

    corrband_mat(windowsize, input_sigma_cor, output_corrband);

    std::cout << "Program End" << std::endl;
    input_x.close();
    input_y.close();
    output_corrband.close();
    input_sigma_cor.close();
    return 0;
}