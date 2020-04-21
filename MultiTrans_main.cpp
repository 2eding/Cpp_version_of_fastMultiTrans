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
    int windowsize = 1000;
    //std::ofstream output_c("result.txt");
    std::vector<double> vc_1;
    std::vector<double> vc_2;

    if (false) {
        std::cout << "This is test verison. please input x, y, windowsize, outputpath" << std::endl;
    }
    //std::cout << "test1" << std::endl;
    //Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
    //std::cout << "test2" << std::endl;
    //Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));
    //std::cout << "test3" << std::endl;
    //Eigen::MatrixXd K = estimateKinship(X);
    //std::cout << "test4" << std::endl;
    //estimateVarComp(K, X, Y, vc_1, vc_2); //calculate vc_1, vc_2
    //std::cout << "test5" << std::endl;
    //Eigen::MatrixXd cor_matrix = cor_mat(K, X, vc_1, vc_2);
    //std::cout << "test6" << std::endl;
    //Eigen::MatrixXd cor_band_mat = corrband_mat(windowsize, cor_matrix);
    ////std::cout << VC <<std::endl;
    //output_c << cor_band_mat << std::endl;
    std::cout << "test7" << std::endl;
    slide_1prep(windowsize, "result.txt", "./test2/prep");
    std::cout << "test8" << std::endl;
    slide_2run("./test2/prep", "./test2/maxstat.txt");
    std::cout << "test9" << std::endl;
    slide_3sort("./test2/maxstat.txt", "./test2/sorted");
    std::cout << "test10" << std::endl;
    slide_4correct("./test2/sorted", "./test2/mul.txt");
    std::cout << "test11" << std::endl;

    vc_1.~vector();
    vc_2.~vector();
    input_x.close();
    input_y.close();
    //output_c.close();

    return 0;
}


