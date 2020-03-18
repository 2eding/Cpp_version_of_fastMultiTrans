#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

void estimateKinship(std::string SNP);
int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);

int main() {
    std::cout << "@----------------------------------------------------------@" << std::endl;
    std::cout << "| MTVCP			|		v1.0		|	10/12/2019		|" << std::endl;
    std::cout << "| (C) 2019 Gi Ju Lee, Tae Gun Kim and Jong Wha Joanne Joo	|" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| For documentation, citation & bug-report instructions:   |" << std::endl;
    std::cout << "|				https://github.com/2eding/MTVCP				|" << std::endl;
    std::cout << "@----------------------------------------------------------@" << std::endl;

    estimateKinship("X.txt");

    return 0;
}
void estimateKinship(std::string SNP) {
    std::ifstream in("X.txt"); // input SNP file
    std::ofstream out("K.txt"); // output kinship file
    std::string read_buffer; // SNP read buffer
    std::string token;
    std::stringstream stream;

    int snp_mat_row = count_matrix_row(in);
    int snp_mat_col = count_matrix_col(in);

    MatrixXd W = MatrixXd(snp_mat_row, snp_mat_col);
    MatrixXd K = MatrixXd(snp_mat_col, snp_mat_col);

    for (int i = 0; i < snp_mat_row; i++) { //row
        std::getline(in, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < snp_mat_col; j++) { //col
            stream >> token;
            W(i, j) = std::stold(token);
        }
        stream.clear();
    }

    W = W.transpose().eval();
    int n = W.rows();
    int m = W.cols();

    // Calculate mean
    RowVectorXd mean = W.colwise().mean();
    MatrixXd matrixMean = mean;

    // Calculate variance
    RowVectorXd var = (W.rowwise() - mean).array().square().colwise().mean();
    MatrixXd matrixVar = var;

    // Calculate standardDeviation
    MatrixXd matrixStd = matrixVar.array().sqrt();
    
    MatrixXd kinship = MatrixXd(m, m);

    // Standardization
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            W(j, i) = (W(j, i) - matrixMean.col(i).value()) / matrixStd.col(i).value();
        }
    }

    // Estimate kinship
    // X * X^T / n
    kinship = W * W.transpose().eval() * 1.0 / m;

    out << kinship << "\n";

    in.close();
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
