#ifndef _INPUT_TEMPLATE_H_
#define _INPUT_TEMPLATE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
using namespace std;
template<typename T>
using matrix = std::vector<std::vector<T>>;

int ReadTridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols);
int ReadMatrixMarketFile(matrix<double>& A, int& num_rows, int& num_cols);
int ReadVectorMarketFile(vector<double>& B, int& num_cols);
int Write1DArrayToMatrixMarketFile(const double B[], int num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols);
int WriteMatrixMarketFile(const matrix<double>& A, int num_rows, int num_cols);
int WriteVectorMarketFile(const std::vector<double>& X, int num_rows, int num_cols);
int WriteTridiagonalToVectorMarketFile(const matrix<double>& A, double B[], vector<double>& X, int num_rows, int num_cols); 
void PrintVector(const std::vector<double>& vector);
void PrintVector2D(const matrix<double>& A); 
std::vector<double> CreateZeroVector(int num_rows);

#endif /* _INPUT_TEMPLATE_H_ */
