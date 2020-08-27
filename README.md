# Cpp_version_of_fastMultiTrans
C++ version of MultiTrans

# Prerequistes for running Source Code
(suggest using Linux environment)

#Basic instructions
```
g++ -std=c++11 -I/home/fastMultiTrans/Eigen MultiTrans_main.cpp MultiTrans_define.cpp -fopenmp -o fastMultiTrans
chmod 755 fastMultiTrans
./fastMultiTrans genotypeData, PhenotypeData, corrbandmatrix_outputpath, correlation_output_path, windowsize
```
