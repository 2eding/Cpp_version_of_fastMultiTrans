#include <iostream>
#include <fstream>
#include "MultiTrans.h"


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
    std::ofstream output_c("result.txt");

    Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
    Eigen::MatrixXd K = estimateKinship(X);

    input_x.close();
    output_c.close();

    return 0;
}


Eigen::MatrixXd estimateKinship(Eigen::MatrixXd snp) {
    std::ofstream out("Ktt.txt"); // output kinship file
    std::string read_buffer; // SNP read buffer
    std::string token;
    std::stringstream stream;



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