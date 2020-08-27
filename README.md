# Cpp_version_of_fastMultiTrans
C++ version of MultiTrans (http://genetics.cs.ucla.edu/multiTrans/)

# Prerequistes for running Source Code
(suggest using Linux environment)

#Basic instructions
```
1. g++ -std=c++11 -I/home/fastMultiTrans/Eigen MultiTrans_main.cpp MultiTrans_define.cpp -fopenmp -o fastMultiTrans

2. chmod 755 fastMultiTrans

3. ./fastMultiTrans genotypeData, PhenotypeData, corrbandmatrix_outputpath, correlation_output_path, windowsize
```
