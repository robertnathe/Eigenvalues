#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <complex>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

template<typename T>
using Matrix2D = vector<vector<T>>;

int writeMatrixMarketMatrix(const MatrixXd& matrix, const string& filename = "A_out.dat");
int writeMatrixMarketVector(const VectorXd& vector, const string& filename = "X_out.dat");
int readMatrixMarketMatrix(MatrixXd& matrix, const string& filename = "A.dat");
int readMatrixMarketVector(VectorXd& vector, const string& filename = "B.dat");
void printExecutionTime(const chrono::time_point<chrono::system_clock>& start, const chrono::time_point<chrono::system_clock>& end);
void printVector(const string& label, const VectorXd& vector);
void printMatrix(const string& label, const MatrixXd& matrix);
int computeArnoldiEigenvalues();
int Arnoldi(const MatrixXd& A, const VectorXd& v0, int m, MatrixXd& V, MatrixXd& H);
int computeDominantEigenvalue();
int inversePowerEigenvalue();
int computeRayleighEigenvalue();
int qrEigenvalue();
double rayleigh_quotient(const Eigen::MatrixXd& A, const Eigen::VectorXd& x);
Eigen::VectorXd rayleigh_quotient_iteration(const Eigen::MatrixXd& A, Eigen::VectorXd x, int max_iter = 100, double tolerance = 1e-6);
int calculateDominantEigenvalue(MatrixXd& matrix, VectorXd& vector, VectorXd& resultVector, double& eigenvalue, double tolerance = 1e-5, size_t maxIterations = 1000);
int eigenvaluePowerMethod(MatrixXd& matrix, double& eigenvalue, VectorXd& resultVector, const VectorXd& initialVector, double tolerance, size_t maxIterations);
void matrixVectorMultiply(const MatrixXd& matrix, const VectorXd& vector, VectorXd& result);
double calculateVectorNorm(const VectorXd& vector);
double calculateInnerProduct(const VectorXd& firstVector, const VectorXd& secondVector);
void scaleVector(VectorXd& vector, double scalar);

int writeMatrixMarketMatrix(const MatrixXd& matrix, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }
    outfile << "%%MatrixMarket matrix coordinate real general\n";
    outfile << matrix.rows() << " " << matrix.cols() << " "; // Number of non-zero entries (worst-case)
    size_t nonZeroCount = 0;
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int col = 0; col < matrix.cols(); ++col) {
            if (matrix(row, col) != 0.0) {
                nonZeroCount++;
            }
        }
    }
    outfile << nonZeroCount << "\n";
    outfile << fixed << setprecision(15);
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int col = 0; col < matrix.cols(); ++col) {
            if (matrix(row, col) != 0.0) {
                outfile << row + 1 << " " << col + 1 << " " << matrix(row, col) << "\n";
            }
        }
    }
    outfile.close();
    return 0;
}

int writeMatrixMarketVector(const VectorXd& vector, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }
    outfile << "%%MatrixMarket matrix coordinate real general\n";
    outfile << 1 << " " << vector.size() << " " << vector.size() << "\n"; // Always a column vector
    outfile << fixed << setprecision(15);
    for (int index = 0; index < vector.size(); ++index) {
        outfile << 1 << " " << index + 1 << " " << vector(index) << "\n";
    }
    outfile.close();
    return 0;
}

int readMatrixMarketMatrix(MatrixXd& matrix, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) return -1;
    string line;
    while (getline(file, line) && line[0] == '%'); // Skip comments
    int numRows, numCols, numEntries;
    istringstream headerStream(line);
    if (!(headerStream >> numRows >> numCols >> numEntries)) return -1;
    matrix.resize(numRows, numCols);
    for (int i = 0; i < numEntries; ++i) {
        int row, col;
        double value;
        if (!(file >> row >> col >> value)) return -1;
        matrix(row - 1, col - 1) = value; // Adjust indices, use Eigen's operator()
    }
    return 0;
}

int readMatrixMarketVector(VectorXd& vector, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening vector file " << filename << endl;
        return -1;
    }
    string line;
    while (getline(file, line) && line[0] == '%');
    int numRows, numCols, numEntries;
    istringstream headerStream(line);
    if (!(headerStream >> numRows >> numCols >> numEntries)) return -1;
    if (numCols != 1 && numRows != 1) {
        cerr << "Error: File does not represent a vector." << endl;
        return -1;
    }
    int vectorSize = (numCols == 1) ? numRows : numCols;
    vector.resize(vectorSize);
    for (int i = 0; i < numEntries; ++i) {
        int row, col;
        double value;
        if (!(file >> row >> col >> value)) return -1;
        vector((numCols == 1) ? row - 1 : col - 1) = value; 
    }
    return 0;
}

void printExecutionTime(const chrono::time_point<chrono::system_clock>& start, const chrono::time_point<chrono::system_clock>& end) {
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << "s\n";
}

void printVector(const string& label, const VectorXd& vector) {
    cout << label << ":\n[";
    cout << fixed << setprecision(5);
    for (int index = 0; index < vector.size(); ++index) {
        cout << vector(index);
        if (index < vector.size() - 1) cout << ", ";
    }
    cout << "]\n" << endl;
}

void printMatrix(const string& label, const MatrixXd& matrix) {
    cout << label << ":\n";
    cout << fixed << setprecision(5);
    for (int row = 0; row < matrix.rows(); ++row) {
        for (int col = 0; col < matrix.cols(); ++col) {
            cout << setw(10) << matrix(row, col);
        }
        cout << endl;
    }
    cout << endl;
}

int computeArnoldiEigenvalues() {
    cout << "Arnoldi Krylov Subspace Algorithm: " << endl << endl;
    MatrixXd A;
    VectorXd v0;
    if (readMatrixMarketMatrix(A, "A.dat") != 0) {
        cerr << "Failed to read matrix A." << endl;
        return -1;
    }
    if (readMatrixMarketVector(v0, "B.dat") != 0) {
        cerr << "Failed to read initial vector." << endl;
        return -1;
    }
    if (A.rows() != v0.size()) {
        cerr << "Error: Matrix A and initial vector dimensions do not match." << endl;
        return -1;
    }
    int m = A.rows();
    MatrixXd V;
    MatrixXd H;
    auto startTime = chrono::system_clock::now();
    int result = Arnoldi(A, v0, m, V, H);
    if (result != 0) {
        cerr << "Arnoldi iteration failed." << endl;
        return result;
    }
    // Extract the upper Hessenberg matrix Hm (m x m)
    MatrixXd Hm = H.topLeftCorner(m, m);
    Eigen::EigenSolver<MatrixXd> eigensolver(Hm);
    if (eigensolver.info() != Eigen::Success) {
        cerr << "Eigenvalue computation failed." << endl;
        return -1;
    }
    VectorXcd eigenvalues = eigensolver.eigenvalues();
    auto endTime = chrono::system_clock::now();
    printMatrix("Hessenberg Matrix H", Hm);
    cout << "Eigenvalues (Ritz values):\n";
    cout << fixed << setprecision(10);
    for (const auto& eig : eigenvalues) {
        cout << eig << endl;
    }
    printExecutionTime(startTime, endTime);
    cout << endl;
    return 0;
}

int Arnoldi(const MatrixXd& A, const VectorXd& v0, int m, MatrixXd& V, MatrixXd& H) {
    int n = A.rows();
    if (n != A.cols()) {
        cerr << "Error: Matrix A must be square." << endl;
        return -1;
    }
    if (v0.size() != n) {
        cerr << "Error: Initial vector has incorrect size." << endl;
        return -1;
    }
    V.resize(n, m + 1);
    H.resize(m + 1, m);
    H.setZero();
    VectorXd v = v0.normalized();
    V.col(0) = v;
    for (int j = 0; j < m; ++j) {
        VectorXd w = A * V.col(j);
        for (int i = 0; i <= j; ++i) {
            H(i, j) = V.col(i).dot(w);
            w -= H(i, j) * V.col(i);
        }
        double h_j1j = w.norm();
        H(j + 1, j) = h_j1j;
        if (h_j1j < 1e-12) {
            if (j < m - 1) {
                cerr << "Breakdown at step " << j + 1 << endl;
                return -1;
            } else {
                H(j + 1, j) = 0.0;
            }
        } else {
            V.col(j + 1) = w / h_j1j;
        }
    }
    return 0;
}

int computeDominantEigenvalue() {
	std::cout << "Power iteration algorithm: " << std::endl << std::endl;
    MatrixXd inputMatrix;
    VectorXd inputVector;
    VectorXd resultVector;
    if (readMatrixMarketMatrix(inputMatrix) != 0) {
        cerr << "Failed to read matrix." << endl;
        return -1;
    }
    if (readMatrixMarketVector(inputVector) != 0) {
        cerr << "Failed to read vector." << endl;
        return -1;
    }
    if (inputMatrix.rows() != inputVector.size()) {
        cerr << "Error: Matrix and vector dimensions do not match." << endl;
        return -1;
    }
    double dominantEigenvalue;
    resultVector.resize(inputMatrix.rows()); // Initialize resultVector
    resultVector.setZero();
    printMatrix("Input Matrix", inputMatrix);
    printVector("Input Vector", inputVector);
    auto startTime = chrono::system_clock::now();
    int result = calculateDominantEigenvalue(inputMatrix, inputVector, resultVector, dominantEigenvalue);
    if (result != 0) {
        cerr << "Error calculating dominant eigenvalue: " << result << endl;
        return result;
    }
    auto endTime = chrono::system_clock::now();
    cout << "Dominant Eigenvalue: " << fixed << setprecision(10) << dominantEigenvalue << endl;
    printVector("Result Vector", resultVector);
    writeMatrixMarketVector(resultVector);
    printExecutionTime(startTime, endTime);
    std::cout << std::endl;
    return 0;
}

int inversePowerEigenvalue() {
	std::cout << "Inverse power iteration algorithm: " << std::endl << std::endl;
    MatrixXd inputMatrix;
    VectorXd inputVector;
    VectorXd resultVector;
    if (readMatrixMarketMatrix(inputMatrix) != 0) {
        cerr << "Failed to read matrix." << endl;
        return -1;
    }
    if (readMatrixMarketVector(inputVector) != 0) {
        cerr << "Failed to read vector." << endl;
        return -1;
    }
    if (inputMatrix.rows() != inputVector.size()) {
        cerr << "Error: Matrix and vector dimensions do not match." << endl;
        return -1;
    }
    double inverseEigenvalue;
    resultVector.resize(inputMatrix.rows());
    resultVector.setZero();
    printMatrix("Input Matrix", inputMatrix);
    printVector("Input Vector", inputVector);
    auto startTime = chrono::system_clock::now();
    MatrixXd pseudoInverseMatrix = inputMatrix.completeOrthogonalDecomposition().pseudoInverse();
    int result = calculateDominantEigenvalue(pseudoInverseMatrix, inputVector, resultVector, inverseEigenvalue);
    if (result != 0) {
        cerr << "Error calculating inverse dominant eigenvalue: " << result << endl;
        return result;
    }
    auto endTime = chrono::system_clock::now();
    if (inverseEigenvalue == 0.0) {
        cerr << "Error: Eigenvalue is zero. Cannot compute inverse." << endl;
        return -1;
    }
    cout << "Inverse Dominant Eigenvalue: " << fixed << setprecision(10) << 1.0 / inverseEigenvalue << endl; // More precise inverse calculation
    printVector("Result Vector", resultVector);
    writeMatrixMarketVector(resultVector);
    printExecutionTime(startTime, endTime);
    std::cout << std::endl;
    return 0;
}

int computeRayleighEigenvalue() {
	std::cout << "Rayleigh iteration algorithm: " << std::endl << std::endl;
    MatrixXd inputMatrix;
    VectorXd initialVector;
    VectorXd intermediateVector;
    if (readMatrixMarketMatrix(inputMatrix) != 0) {
        cerr << "Failed to read input matrix." << endl;
        return -1;
    }
    if (readMatrixMarketVector(initialVector) != 0) {
        cerr << "Failed to read initial vector." << endl;
        return -1;
    }
    if (inputMatrix.rows() != initialVector.size()) {
        cerr << "Error: Input matrix and vector dimensions do not match." << endl;
        return -1;
    }
    double dominantRayleighEigenvalue;
    intermediateVector.resize(inputMatrix.rows());
    intermediateVector.setZero();
    printMatrix("Input Matrix", inputMatrix);
    printVector("Initial Vector", initialVector);
    auto startTime = chrono::system_clock::now();
    initialVector = rayleigh_quotient_iteration(inputMatrix, initialVector, 1000, 1e-8);
    dominantRayleighEigenvalue = rayleigh_quotient(inputMatrix, initialVector);
    cout << "Dominant Rayleigh Eigenvalue: " << fixed << setprecision(10) << dominantRayleighEigenvalue << endl;
    printVector("Updated Vector", initialVector);
    writeMatrixMarketVector(initialVector);
    auto endTime = chrono::system_clock::now();
    printExecutionTime(startTime, endTime);
    std::cout << std::endl;
    return 0;
}

int qrEigenvalue() {
	std::cout << "QR Hessenberg direct algorithm: " << std::endl << std::endl;
    MatrixXd inputMatrix;
    if (readMatrixMarketMatrix(inputMatrix) != 0) {
        cerr << "Failed to read matrix." << endl;
        return -1;
    }
    printMatrix("Input Matrix", inputMatrix);
    auto startTime = chrono::system_clock::now();
    // 1. Direct QR decomposition using Eigen
    Eigen::ComplexEigenSolver<MatrixXd> eigenSolver(inputMatrix); // Use ComplexEigenSolver for general case
    if (eigenSolver.info() != Eigen::Success) {
        cerr << "Eigenvalue decomposition failed!" << endl;
        return -1;
    }
    VectorXcd eigenvalues = eigenSolver.eigenvalues();
    MatrixXcd eigenvectors = eigenSolver.eigenvectors();
    auto endTime = chrono::system_clock::now();
    // 2. Print results (more robust printing for complex numbers)
    cout << "Eigenvalues:\n";
    cout << fixed << setprecision(10);
    for (int index = 0; index < eigenvalues.size(); ++index) {
        cout << eigenvalues(index) << "\n";
    }
    cout << endl;
    cout << "Eigenvectors:\n";
    for (int index = 0; index < eigenvectors.cols(); ++index) {
        cout << "Eigenvector " << index + 1 << ":\n";
        for (int row = 0; row < eigenvectors.rows(); ++row) {
            cout << eigenvectors(row, index) << "\n";
        }
        cout << endl;
    }
    printExecutionTime(startTime, endTime); 
    return 0;
}

double rayleigh_quotient(const Eigen::MatrixXd& A, const Eigen::VectorXd& x) {
  double numerator = x.transpose() * A * x;
  double denominator = x.transpose() * x;
  return numerator / denominator;
}

Eigen::VectorXd rayleigh_quotient_iteration(const Eigen::MatrixXd& A, Eigen::VectorXd x, int max_iter, double tolerance) {
  double prev_eigenvalue = 0.0;
  for (int i = 0; i < max_iter; ++i) {
    double eigenvalue = rayleigh_quotient(A, x);
    if (std::abs(eigenvalue - prev_eigenvalue) < tolerance) {
      break;
    }
    prev_eigenvalue = eigenvalue;
    Eigen::VectorXd w = (A - eigenvalue * Eigen::MatrixXd::Identity(A.rows(), A.cols())).lu().solve(x);
    x = w / w.norm();
  }
  return x;
}

int calculateDominantEigenvalue(MatrixXd& matrix, VectorXd& vector, VectorXd& resultVector, double& eigenvalue, double tolerance, size_t maxIterations) {
    return eigenvaluePowerMethod(matrix, eigenvalue, resultVector, vector, tolerance, maxIterations);
}

int eigenvaluePowerMethod(MatrixXd& matrix, double& eigenvalue, VectorXd& resultVector, const VectorXd& initialVector, double tolerance, size_t maxIterations) {
    double previousEigenvalue = 0.0;
    resultVector = initialVector;
    for (size_t iteration = 0; iteration < maxIterations; ++iteration) {
        matrixVectorMultiply(matrix, resultVector, resultVector); // resultVector = matrix * resultVector (in-place)
        eigenvalue = calculateVectorNorm(resultVector);
        if (eigenvalue == 0.0) {
            return -2; // Iteration vector vanished
        }
        scaleVector(resultVector, eigenvalue); // Normalize resultVector
        if (iteration > 0 && abs((eigenvalue - previousEigenvalue) / eigenvalue) < tolerance) {
            return 0; // Convergence achieved
        }
        previousEigenvalue = eigenvalue;
    }
    return -1; // Iteration did not converge
}

void matrixVectorMultiply(const MatrixXd& matrix, const VectorXd& vector, VectorXd& result) {
    result.resize(matrix.rows()); // Make sure result is the correct size
    result = matrix * vector;          // Use Eigen's built-in multiplication
}

double calculateVectorNorm(const VectorXd& vector) {
    return vector.norm();    // Use Eigen's built-in norm calculation
}

double calculateInnerProduct(const VectorXd& firstVector, const VectorXd& secondVector) {
    return firstVector.dot(secondVector);    // Use Eigen's built-in dot product
}

void scaleVector(VectorXd& vector, double scalar) {
    vector /= scalar;              // Use Eigen's built-in scaling
}

int main() {
    computeDominantEigenvalue(); 
    inversePowerEigenvalue();
    computeRayleighEigenvalue();
    qrEigenvalue();
    computeArnoldiEigenvalues();
    return 0;
}
