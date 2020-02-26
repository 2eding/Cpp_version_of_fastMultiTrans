
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Eigenvalues"
#include "Eigen/unsupported/Eigen/MatrixFunctions"



void generateCovBand(int windowSize, std::string corrMatrix, std::string output);
void generateCovR(std::string X_rightdim, std::string Kpath, std::string VCpath, std::string outputPath);
int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);
double cal_median(std::vector<std::string> col);
void cal_standard_deviation(std::vector<double>& ret, Eigen::MatrixXd& mat);
void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat, Eigen::MatrixXd& colmat);

//int argc, char *argv[]

int main(int argc, char* argv[])
{
    //if (argc == 4) {
    //   generateCovBand(atoi(argv[1]), argv[2], argv[3]);
    //}
    //else {
    //   std::cout << "Usage: generateCovBand [window_size] [inputPath] [outputPath]" << std::endl;
    //}
    //return 0;

    generateCovR("X_rightdim.txt", "K.txt", "VC.txt", "test.txt");

    return 0;
}

void generateCovR(std::string X_rightdim, std::string Kpath, std::string VCpath, std::string outputPath) {
    std::ifstream in_X_rightdim(X_rightdim.c_str()); //input X_rightdim matrix file
    std::ifstream in_Kpath(Kpath.c_str()); //input Kpath matrix file
    std::ifstream in_VCpath(VCpath.c_str()); //input VCpath matrix file
    std::ofstream out(outputPath); //output path

    std::string read_buffer; //file preprocessing variable
    std::string token;
    std::stringstream stream;

    Eigen::MatrixXd VC_matrix;
    

    std::vector<std::string> VC_col1;
    std::vector<std::string> VC_col2;

    std::vector<double> stand_devi;

    Eigen::MatrixXd X_right_matrix;

    Eigen::MatrixXd K_matrix;
    Eigen::MatrixXd K_matrix_eigen_vectors;
    Eigen::MatrixXd K_matrix_eigen_values;

    Eigen::MatrixXd U_matrix;
    Eigen::MatrixXd cor_matrix;

    Eigen::MatrixXd test = Eigen::MatrixXd(2, 2);

    /*for (int i = 0; i < 3; i++) {
       for (int j = 0; j < 3; j++) {
          test(j, i) = i;
       }
    }*/
    //test(0, 0) = 1; test(0, 1) = 3; test(0, 2) = 7; 
    //test(1, 0) = 2; test(1, 1) = 4; test(1, 2) = 8; 
   // test(2, 0) = 3; test(2, 1) = 6; test(2, 2) = 9; 
    //Eigen::EigenSolver<Eigen::MatrixXd> tt(test);
    //std::cout << test << std::endl;
    //std::cout << tt.pseudoEigenvectors() << std::endl << std::endl;
    //std::cout << tt.pseudoEigenvalueMatrix() << std::endl << std::endl;

    //std::cout << test * tt.pseudoEigenvectors().col(0) << std::endl << std::endl;
    //std::cout << tt.pseudoEigenvalueMatrix()(0,0) * tt.pseudoEigenvectors().col(0) << std::endl;
    //
    //std::cout << std::endl;
    //std::cout << test * tt.pseudoEigenvectors().col(1) << std::endl << std::endl;
    //std::cout << tt.pseudoEigenvalueMatrix()(1, 1) * tt.pseudoEigenvectors().col(1) << std::endl;

    //std::cout << std::endl;
    //std::cout << test * tt.pseudoEigenvectors().col(2) << std::endl << std::endl;
    //std::cout << tt.pseudoEigenvalueMatrix()(2, 2) * tt.pseudoEigenvectors().col(2) << std::endl;
    //std::cout << test * test << std::endl;
    //exit(1);

    int snpNum = 0; int indiNum = 0; //X_rightdim's row, col number
    double Median_Vg = 0; double Median_Ve = 0; //Vc's Median like median(VC[,1]), median(VC[,2]

    snpNum = count_matrix_col(in_X_rightdim); //check file matrix size
    indiNum = count_matrix_row(in_X_rightdim);
    int k_size = count_matrix_col(in_Kpath);

    //std::cout << snpNum << std::endl;
    //std::cout << indiNum << std::endl;

    if (snpNum == 0 || indiNum == 0) {
        std::cout << "X_rightdim is empty, please check input file or path" << std::endl;
        exit(1);
    }
    else if (k_size != count_matrix_row(in_Kpath)) {
        std::cout << "K file is not n by n matrix" << std::endl;
        exit(1);
    }

    while (in_VCpath.peek() != EOF) {
        std::getline(in_VCpath, read_buffer);
        stream.str(read_buffer);
        stream >> token;  //std::cout << "token1 : " << token <<std::endl;
        VC_col1.push_back(token);
        stream >> token; //std::cout << "token2 : " << token <<std::endl; 
        VC_col2.push_back(token);
        stream.clear();
    }


    Median_Vg = cal_median(VC_col1); Median_Ve = cal_median(VC_col2);  //cal median value
    VC_col1.clear(); VC_col2.clear();

    std::cout << "Median_Vg = " << Median_Vg << std::endl;
    std::cout << "Median_Ve = " << Median_Ve << std::endl;

    X_right_matrix = Eigen::MatrixXd(indiNum, snpNum);
    K_matrix = Eigen::MatrixXd(k_size, k_size);
    K_matrix_eigen_vectors = Eigen::MatrixXd(k_size, k_size);
    K_matrix_eigen_values = Eigen::MatrixXd(k_size, k_size);
    U_matrix = Eigen::MatrixXd(k_size, k_size);
    cor_matrix = Eigen::MatrixXd(snpNum, snpNum);
    //make X_matrix
    for (int i = 0; i < indiNum; i++) { //row
        std::getline(in_X_rightdim, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < snpNum; j++) { //col
            stream >> token;
            X_right_matrix(i, j) = std::stold(token);

            //std::cout.precision(16); 
            //std::cout << X_right_matrix(i,j) << " ";
        }
        stream.clear();
        //std::cout << "" << std::endl;
    }

    // make K_matrix
    for (int i = 0; i < k_size; i++) { //row
        std::getline(in_Kpath, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < k_size; j++) { //col
            stream >> token;
            if (i != j) {
                K_matrix(i, j) = Median_Vg * std::stold(token);
            }
            else { //i == j
                K_matrix(i, j) = (Median_Vg * std::stold(token)) + Median_Ve;
            }
            //std::cout.precision(16); // check 
            //std::cout << K_matrix(i,j) << " ";
        }
        //stream.clear();
        //std::cout << "" << std::endl;
    } //then K_matrix -> sigmaM matrix
    //exit(1);
    Eigen::EigenSolver<Eigen::MatrixXd> cal_eigen(K_matrix);
    K_matrix_eigen_vectors = cal_eigen.pseudoEigenvectors();
    K_matrix_eigen_values = cal_eigen.pseudoEigenvalueMatrix();
    
    K_matrix = K_matrix.sqrt();
    //for (int i = 0; i < K_matrix_eigen_vectors.rows(); i++) {
    //    for (int j = 0; j < K_matrix_eigen_vectors.cols(); j++) {
    //        if (j + 1 == K_matrix_eigen_vectors.cols()) {
    //            out << K_matrix_eigen_vectors(i, j);
    //        }
    //        else {
    //            out << K_matrix_eigen_vectors(i, j) << " ";
    //        }
    //    }out << "\n";
    //}

    //for (int i = 0; i < k_size; i++) {
    //        if (K_matrix_eigen_values(i, i) < 0.0000000000001) { // 1e-1
    //            K_matrix_eigen_values(i, i) = 0.0000000000001;
    //        }
    //        K_matrix_eigen_values(i, i) = 1/sqrt(K_matrix_eigen_values(i, i));
    //        //std::cout << i + 1 << ": " << K_matrix_eigen_values(i, i) << std::endl;
    //}
    //U_matrix = K_matrix_eigen_vectors * K_matrix_eigen_values;
    //
    ////std::cout << K_matrix_eigen_vectors << std::endl << std::endl;
    //K_matrix_eigen_vectors.transposeInPlace();
    //U_matrix = U_matrix * K_matrix_eigen_vectors;
    //U_matrix.transposeInPlace();
    //std::cout << U_matrix << std::endl;
    K_matrix.transposeInPlace();
    K_matrix = K_matrix * X_right_matrix;
    cal_standard_deviation(stand_devi, K_matrix);
    cal_cor(stand_devi, K_matrix, cor_matrix);
    //U_matrix = U_matrix * X_right_matrix;
    //cal_standard_deviation(stand_devi, U_matrix);
    //std::cout << U_matrix << std::endl;
    //for (int i = 0; i < U_matrix.rows(); i++) {
    //    for (int j = 0; j < U_matrix.cols(); j++) {
    //        if (j + 1 == U_matrix.cols()) {
    //            out << U_matrix(i, j);
    //        }
    //        else {
    //            out << U_matrix(i, j) << " ";
    //        }
    //    }out << "\n";
    //}
    //cal_cor(stand_devi, U_matrix, cor_matrix);
    
    out.precision(16);
    out << cor_matrix << std::endl;

    stand_devi.~vector();
    
    cor_matrix.resize(0, 0);
    X_right_matrix.resize(0, 0);
    K_matrix_eigen_values.resize(0, 0);
    K_matrix_eigen_vectors.resize(0, 0);
    U_matrix.resize(0, 0);

    in_X_rightdim.close();
    in_Kpath.close();
    in_VCpath.close();
    out.close();
}

void generateCovBand(int windowSize, std::string corrMatrix, std::string output) {
    std::ifstream in(corrMatrix.c_str()); //input matrix file
    std::ofstream out(output); //output path
    std::string read_buffer; //corrMatrix read buffer
    std::getline(in, read_buffer); //read first line(perhaps.. preprocessing)
    int snpCnt = 0; //snip count??
    std::string token;
    std::stringstream stream;

    while (in.peek() != EOF) {
        ++snpCnt;
        std::getline(in, read_buffer);
stream.str(read_buffer);

for (int i = 0; i < windowSize - snpCnt; i++) {
    out << "0\t";
}

for (int i = 0; i < snpCnt - windowSize; i++) {
    stream >> token;
}

if (snpCnt <= windowSize) {
    for (int i = 0; i < snpCnt; i++) {
        stream >> token;
        out << token << "\t";
    }
}
else {
    for (int i = 0; i < windowSize; i++) {
        stream >> token;
        out << token << "\t";
    }
}
out << "\n";
    } //end while
    std::cout << snpCnt << " number of lines read" << std::endl;
    in.close();
    out.close();

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

double cal_median(std::vector<std::string> col) {
    std::vector<double> temp;
    double median = 0.0;
    for (int i = 0; i < col.size(); i++) {
        temp.push_back(std::stod(col[i]));
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
    double temp;
    double temp_sum;
    for (int i = 0; i < mat.cols(); i++) { //col
        temp = 0;
        temp_sum = 0;
        for (int j = 0; j < mat.rows(); j++) { //row
            temp = mat(j, i) - (mat.col(i).sum() / ((double)mat.rows() - 1.0));
            temp_sum += temp*temp; // mat.col(i).sum()/(double)mat.rows() is col average
        }
        //std::cout << sqrt((temp_sum) / ((double)mat.rows() - 1.0)) << std::endl;
        ret.push_back(sqrt((temp_sum)));

    }
}

void cal_cor(std::vector<double>& ret, Eigen::MatrixXd& mat, Eigen::MatrixXd& colmat) {
    double temp;
    double total_temp;
    for (int k = 0; k < mat.cols(); k++) { // cols
        for (int i = 0; i < mat.cols(); i++) { // cols
            total_temp = 0;
            for (int j = 0; j < mat.rows(); j++) { // rows 
                temp = mat(j, k) - (mat.col(k).sum() / (double)(mat.rows() - 1.0));
                temp = temp * (mat(j, i) - (mat.col(i).sum() / (double)(mat.rows() - 1.0)));
/*                std::cout << "col "<< i <<" avr :" << (mat.col(i).sum() / (double)(mat.rows() - 1.0)) << std::endl;
                std::cout << "col "<< k <<" avr :" << (mat.col(k).sum() / (double)(mat.rows() - 1.0)) << std::endl;
                std::cout << j << " :" << temp << std::endl;
           */     total_temp += temp;
            }
            colmat(k, i) = (((total_temp))/(ret[k]*ret[i]));
        }
    }
}
