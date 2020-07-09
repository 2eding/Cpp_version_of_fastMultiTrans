#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/Core>

void estimateKinship(std::string SNP);
Eigen::MatrixXd calculateKinship(Eigen::MatrixXd W);
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

    estimateKinship("X2.txt");

    return 0;
}

void estimateKinship(std::string SNP) {
    std::ifstream in("X2.txt"); // input SNP file
    std::ofstream out("K22.txt"); // output kinship file
    std::string read_buffer; // SNP read buffer
    std::string token;
    std::stringstream stream;

    std::clock_t start = std::clock();

    int snp_mat_row = count_matrix_row(in);
    int snp_mat_col = count_matrix_col(in);
    Eigen::MatrixXd X = Eigen::MatrixXd(snp_mat_row, snp_mat_col);

    int n = snp_mat_col;
    int m = 1000;
    Eigen::MatrixXd W = Eigen::MatrixXd(n, m).setOnes() * NAN;

    int i = 0;
    Eigen::MatrixXd kinship = Eigen::MatrixXd(n, n).setZero();
    Eigen::MatrixXd kinship_j = Eigen::MatrixXd(n, n);
    while (i < snp_mat_row) {
        int j = 0;
        while ((j < m) && (i < snp_mat_row)) {
    
            std::getline(in, read_buffer);
            stream.str(read_buffer);
            for (int k = 0; k < snp_mat_col; k++) {
                stream >> token;
                X(i, k) = std::stold(token);
            }

            double mean = X.row(i).mean();
            double snp_var = X.row(i).squaredNorm() - (mean * mean);
            stream.clear();
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

        std::cout << "Processing first " << i << " SNPs" << std::endl;
        
        if (kinship.isZero()) {
            kinship = calculateKinship(W.transpose()) * j;
        }
        else {
            kinship_j = calculateKinship(W.transpose()) * j;
            kinship = kinship + kinship_j;
        }
    }
    kinship = kinship / float(snp_mat_row);
    double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "running time = " << duration << " secs\n";
    out << kinship << "\n";

    in.close();
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
