#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <float.h>                // required for DBL_EPSILON
#include <chrono>
using namespace std;
template<typename T>
using matrix = std::vector<std::vector<T>>;

void PrintMatrix(const vector<double> A[], unsigned int num_rows);
void PrintVector(const double B[], unsigned int num_rows);
void Identity_Matrix(double *U, unsigned int n);
int ReadMatrixMarketMatrix(matrix<double>& A, int& num_rows, int& num_cols);
int ReadMatrixMarketVector(vector<double>& B, int& num_cols);
int Write1DArrayToMatrixMarketFile(const double B[], int num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols);

int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols);
int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols);
void PrintVector(const std::vector<double>& vector);
void PrintMatrix(const matrix<double>& A); 
std::vector<double> InitializeVectorToZero(int num_rows);

/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
void PrintMatrix(const vector<double> A[], unsigned int num_rows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < num_rows; i++) {
    for (
      auto it = A[i].begin();
        it != A[i].end(); it++)
        cout << *it << " ";
      cout << endl;
    }  
}

void PrintVector(const double TempVector[], unsigned int num_rows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < num_rows; i++)
    std::cout << TempVector[i] << "   "; 
  std::cout << endl;
}

////////////////////////////////////////////////////////////////////////////////
// File: identity_matrix.c                                                    //
// Routine(s):                                                                //
//    Identity_Matrix                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Identity_Matrix(double *A, int n)                                    //
//                                                                            //
//  Description:                                                              //
//     Set the square n x n matrix A equal to the identity matrix, i.e.       //
//     A[i][j] = 0 if i != j and A[i][i] = 1.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    n     The number of rows and columns of the matrix A.           //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//                                                                            //
//     Identity_Matrix(&A[0][0], N);                                          //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Identity_Matrix(double *U, unsigned int n)
{
   unsigned int i,j;
   for (i = 0; i < n - 1; i++) {
      *U++ = 1.0;
      for (j = 0; j < n; j++) *U++ = 0.0;
   } 
   *U = 1.0;
}

int ReadMatrixMarketMatrix(matrix<double>& A, int& num_rows, int& num_cols) {

 int number_of_entries_A {0};
 int i_index {0}, j_index {0};
 double elem {0.0};
 FILE* myFile = fopen("A.dat", "r");
 if (myFile == NULL) {
   std::cerr << "Error Reading File" << endl;
   exit(0);
 }
 // Skip header comments
 fscanf(myFile, "%*[^\n]\n"); // Read and discard header line
 // Read dimensions
 if (fscanf(myFile, "%d %d %d\n", &num_rows, &num_cols, &number_of_entries_A) != 3) {
   std::cerr << "Error reading matrix dimensions from A.dat" << endl;
   fclose(myFile);
   return -1;
 }
 // Resize A to accommodate num_num_rows and num_num_num_cols
 A.resize(num_rows);
 for (int i = 0; i < num_rows; ++i) {
   A[i].resize(num_cols);
 }
 // Read non-zero elements by row and column indices
 for (int i = 0; i < number_of_entries_A; ++i) {
   if (fscanf(myFile, "%d %d %lf\n", &i_index, &j_index, &elem) != 3) {
     std::cerr << "Error reading matrix entries from A.dat" << endl;
     fclose(myFile);
     return -1;
   }
   i_index--; // Adjust for zero-based indexing
   j_index--;
   A[i_index][j_index] = elem;
 }
 fclose(myFile);
 return 0;
}

// ReadMatrixMarketVector implementation for vectors
int ReadMatrixMarketVector(vector<double>& B, int& num_cols) {
    (void) num_cols; 
   FILE *myFile2;
   myFile2 = fopen ("B.dat", "r");
   int dim_B_Array[3];
   int i_index {0}, j_index {0};
   double value {0.0};
   while (myFile2 == NULL)
   {
    std::cout << "Error Reading File" << endl;
     exit (0);
   } 
   fscanf (myFile2, "%*s %*s %*s %*s %*s");
   for (int i = 0; i < 3; i++)
   {
     fscanf (myFile2, "%d,", &dim_B_Array[i]);
   }
   for (int i = 0; i < dim_B_Array[1]; i++)
     B.push_back(0.0);
   for (int i = 0; i < dim_B_Array[1]; i++)
   {
     fscanf (myFile2, "%d,", &i_index);
     i_index--;
     fscanf (myFile2, "%d,", &j_index);
     j_index--;
     fscanf (myFile2, "%lf,", &value);
     if (value != 0.0) 
     {
       B[i] = value;
     }
   }
   fclose (myFile2);
 return 0;
}

int Write1DArrayToMatrixMarketFile(const double B[], int num_rows) {  
  //Use C++ streams for safer file handling
  ofstream outfile("B_out.dat");
  if (!outfile.is_open()) {
    cerr << "Error opening file for writing: " << "B_out.dat" << endl;
    return 1;
  }
  
  // Use C++ streams for safer file handling
  //ofstream outfile(filename);
  //if (!outfile.is_open()) {
    //cerr << "Error opening file for writing: " << filename << endl;
    //return 1;
  //}
  
  printf ("B =\n");
  // Write header information (assuming general coordinate pattern)
  outfile << "%%MatrixMarket_Output_vector_B.dat matrix coordinate pattern general\n";
  outfile << num_rows << " 1 " << num_rows << endl; // Adjust for 1D array
  // Write each element with row and column indices (starting from 1)
  for (int i = 0; i < num_rows; i++) {
    outfile << i + 1 << " " << 1 << " " << B[i] << endl;
    printf ("%6.5lf    ", B[i]);
  }
  std::cout << std::endl;
  outfile.close();
  return 0;
}

int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols) {
  // Use C++ streams for safer file handling
  ofstream outfile("A_out.dat");
  if (!outfile.is_open()) {
    cerr << "Error opening file for writing: A_out.dat" << endl;
    return 1;
  }

  // Write header information (assuming general coordinate pattern)
  outfile << "%%MatrixMarket_Output_vector_A.dat matrix coordinate pattern general\n";
  outfile << num_rows << " " << num_cols << " " << num_rows * num_cols << endl;
  printf ("A =\n");
  // Write each element with row and column indices (starting from 1)
  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      outfile << i + 1 << " " << j + 1 << " " << array_A[i * num_cols + j] << endl;
      printf ("%6.5lf    ", array_A[i * num_cols + j]);
    }
    std::cout << std::endl;
  }

  outfile.close();
  return 0;
}

int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols) {
 // Open file in writing mode
 int CountNonZeroEntries {0};
 FILE* fptr = fopen("A_out.txt", "w");
 if (fptr == NULL) {
   std::cerr << "Error opening file for writing matrix A." << std::endl;
   return -1; // Indicate error
 }
 for (int i = 0; i < num_rows; ++i) {
   for (int j = 0; j < num_cols; ++j) {
    if (A[i][j] != 0.0)
    {
      CountNonZeroEntries++;
    }
   }
 }
 // Write Matrix Market header information for a sparse matrix
 fprintf(fptr, "%%MatrixMarket_Input_matrix_A.dat matrix coordinate pattern general\n");
// fprintf(fptr, "%d %d %d\n", num_rows, num_cols, num_rows*num_cols); // Count all entries
 fprintf(fptr, "%d %d %d\n", num_rows, num_cols, CountNonZeroEntries); // Count non-zero entries
 // Write only non-zero elements of A
 // I changed this section of code.
 for (int i = 0; i < num_rows; ++i) {
   for (int j = 0; j < num_cols; ++j) {
     if (A[i][j] != 0.0) 
     { // Check for non-zero value
       fprintf(fptr, "%d %d %lf\n", i + 1, j + 1, A[i][j]);
     }
   }
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols) {
 // Open file in writing mode
 FILE* fptr = fopen("X.dat", "w");
 if (fptr == NULL) {
   std::cerr << "Error opening file for writing X vector." << std::endl;
   return -1; // Indicate error
 }
 // Write Matrix Market header information
 fprintf(fptr, "%%MatrixMarket_Input_matrix_X.dat matrix coordinate pattern general\n");
 fprintf(fptr, "%d %d %d\n", 1, num_cols, num_cols); // All entries are assumed non-zero
 std::cout << "%%MatrixMarket_Input_matrix_X.dat matrix coordinate pattern general\n";
 // Write each element of X
 for (int i = 0; i < num_rows; ++i) {
   if (X[i] != 0.0) 
   { // Check for non-zero value
     fprintf(fptr, "%d %d %lf\n", 1, i+1, X[i]); // Row index always 1 for a vector
     std::cout << 1 << "   " << i+1 << "   "<< X[i] << std::endl;
   }
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

void PrintVector(const std::vector<double>& vector) {
// std::cout << "Displaying the vector: " << endl;
 std::cout << std::setprecision(5);
 for (const double& value : vector) {
   std::cout << value << " ";
 }
 std::cout << std::endl;
}

void PrintMatrix(const matrix<double>& A) {
// std::cout << "Displaying the 2D vector:" << endl;
 for (const auto& row : A) {
   for (const double& value : row) {
     std::cout << value << " ";
   }
   std::cout << std::endl;
 }
}

std::vector<double> InitializeVectorToZero(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
}

int QR_Input(double array_A[], double eigen_real[], double eigen_imag[], unsigned int num_rows, unsigned int num_cols, unsigned int Max_Iter);
int Hessenberg_Form_Orthogonal(double *A, double *U, int n);
int QR_Hessenberg_Matrix( double *H, double *S, double eigen_real[],
                          double eigen_imag[], int n, unsigned int Max_Iter);
////////////////////////////////////////////////////////////////////////////////
// File: qr_hessenberg_matrix.c                                               //
// Routine(s):                                                                //
//    QR_Hessenberg_Matrix                                                    //
////////////////////////////////////////////////////////////////////////////////
//                        Internally Defined Routines 
static void   One_Real_Eigenvalue(double Hrow[], double eigen_real[],
                                  double eigen_imag[], int row, double shift);
static void   Two_Eigenvalues(double *H, double *S, double eigen_real[],
                                 double eigen_imag[], int n, int k, double t);
static void   Update_Row(double *Hrow, double cos, double sin, int n, int k);
static void   Update_Column(double* H, double cos, double sin, int n, int k);
static void   Update_Transformation(double *S, double cos, double sin,
                                                               int n, int k);
static void   Double_QR_Iteration(double *H, double *S, int row, int min_row,
                                         int n, double* shift, unsigned int interation);
static void   Product_and_Sum_of_Shifts(double *H, int n, int max_row,
                    double* shift, double *trace, double *det, unsigned int interation);
static int    Two_Consecutive_Small_Subdiagonal(double* H, int min_row,
                                int max_row, int n, double trace, double det);
static void   Double_QR_Step(double *H, int min_row, int max_row, int min_col,
                                  double trace, double det, double *S, int n);
static void   Complex_Division(double x, double y, double u, double v,
                                                       double* a, double* b);
static void   BackSubstitution(double *H, double eigen_real[],
                                                  double eigen_imag[], int n);
static void   BackSubstitute_Real_Vector(double *H, double eigen_real[],
                 double eigen_imag[], int row,  double zero_tolerance, int n);
static void   BackSubstitute_Complex_Vector(double *H, double eigen_real[],
                 double eigen_imag[], int row,  double zero_tolerance, int n);
static void   Calculate_Eigenvectors(double *H, double *S, double eigen_real[],
                                                  double eigen_imag[], int n);
//                    Required Externally Defined Routines 

int QR_Input(double array_A[], double eigen_real[], double eigen_imag[], unsigned int num_rows, unsigned int num_cols, unsigned int Max_Iter)
{
  unsigned int i {0}, j {0};
  int n = num_rows, k {0};
  (void) k;
  (void) n;
  double H[num_rows][num_cols] {0.0};
  double S[num_rows][num_cols] {0.0};
  double U[num_rows][num_cols] {0.0};
  for (i = 0; i < num_rows; i++)
  {  
    for (j = 0; j < num_cols; j++)
	{
	  H[i][j] = array_A[i*num_rows+j];
	  if (i != j){
	    S[i][j] = 0.0;
	    U[i][j] = 0.0;
	  }
	  else {
	    S[i][j] = 1.0;
	    U[i][j] = 1.0;
	  }
	}
  }
  cout << endl;
  Hessenberg_Form_Orthogonal(&H[0][0], *U, n);
  cout << "The Householder Form Orthogonal matrix H is the following: " << endl;
  for (i = 0; i < num_rows; i++)
  {
    for (j = 0; j < num_cols; j++)
	{
//	  cout << S[i][j] << " ";  
	  cout << H[i][j] << " ";  
    }
    cout << endl;
  }
  cout << endl;
  k = QR_Hessenberg_Matrix(*H, *S, &eigen_real[0], &eigen_imag[0], n, Max_Iter);
  if (k < 0) 
  {
	  printf("Failed"); 
	  exit(1);
  }
  cout << "The column eigenvector S matrix is the following: " << endl;
//  cout << "The upper Hessenberg matrix H is the following: " << endl;
  for (i = 0; i < num_rows; i++)
  {
    for (j = 0; j < num_cols; j++)
	{
	  cout << S[i][j] << " ";  
//	  cout << H[i][j] << " ";  
    }
    cout << endl;
  }
  cout << endl;
  cout << "The real part of the eigenvalues are the following: " << endl; 
  PrintVector(&eigen_real[0], num_rows);
  cout << "The imaginary part of the eigenvalues are the following: " << endl;
  PrintVector(&eigen_imag[0], num_rows);
  return 0;	
}

////////////////////////////////////////////////////////////////////////////////
// File: hessenberg_orthog.c                                                  //
// Routine(s):                                                                //
//    Hessenberg_Form_Orthogonal                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Hessenberg_Form_Orthogonal(double *A, double *U, int n)               //
//                                                                            //
//  Description:                                                              //
//     This program transforms the square matrix A to a similar matrix in     //
//     Hessenberg form by a multiplying A on the right and left by a sequence //
//     of Householder transformations.                                        //
//     Def:  Two matrices A and B are said to be orthogonally similar if there//
//           exists an orthogonal matrix U such that A U = U B.               //
//     Def   A Hessenberg matrix is the sum of an upper triangular matrix and //
//           a matrix all of whose components are 0 except possibly on its    //
//           subdiagonal.  A Hessenberg matrix is sometimes said to be almost //
//           upper triangular.                                                //
//     Def:  A Householder transformation is an orthogonal transformation of  //
//           the form Q = I - 2 uu'/u'u, where u is a n x 1 column matrix and //
//           ' denotes the transpose.                                         //
//     Thm:  If Q is a Householder transformation then Q' = Q  and  Q' Q = I, //
//           i.e. Q is a symmetric involution.                                //
//     The algorithm proceeds by successivly selecting columns j = 0,...,n-3  //
//     and then calculating the Householder transformation Q which annihilates//
//     the components below the subdiagonal for that column and leaves the    //
//     previously selected columns invariant.  The algorithm then updates     //
//     the matrix A, in place, by premultiplication by Q followed by          //
//     postmultiplication by Q.                                               //
//     If the j-th column of A is (a[0],...,a[n-1]), then  choose u' =        //
//     (u[0],...,u[n-1]) by u[0] = 0, ... , u[j] = 0, u[j+2] = a[j+2],...,    //
//     u[n-1] = a[n-1].  The remaining component u[j+1] = a[j+1] - s, where   //
//     s^2 = a[j+1]^2 + ... + a[n-1]^2, and the choice of sign for s,         //
//     sign(s) = -sign(a[j+1]) maximizes the number of significant bits for   //
//     u[j+1].                                                                //
//                                                                            //
//     Remark:  If H = U' A U, where U is orthogonal, and if v is an eigen-   //
//     vector of H with eigenvalue x, then A U v = U H v = x U v, so that     //
//     U v is an eigenvector of A with corresponding eigenvalue x.            //
//                                                                            //
//  Arguments:                                                                //
//     double *A   Pointer to the first element of the matrix A[n][n].        //
//                 The original matrix A is replaced with the orthogonally    //
//                 similar matrix in Hessenberg form.                         //
//     double *U   Pointer to the first element of the matrix U[n][n].  The   //
//                 orthogonal matrix which transforms the input matrix to     //
//                 an orthogonally similar matrix in Hessenberg form.         //
//     int     n   The number of rows or columns of the matrix A.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - Unable to allocate space for working storage.             //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], U[N][N];                                               //
//                                                                            //
//     (your code to create the matrix A)                                     //
//     Hessenberg_Form_Orthogonal(&A[0][0], (double*) U, N);                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Hessenberg_Form_Orthogonal(double *A, double *U, int n)
{
   int i, k, col;
   double *u;
   double *p_row, *psubdiag;
   double *pA, *pU;
   double sss;                             // signed sqrt of sum of squares
   double scale;
   double innerproduct;
         // n x n matrices for which n <= 2 are already in Hessenberg form
   Identity_Matrix(U, n);
   if (n <= 2) return 0;
         // Reserve auxillary storage, if unavailable, return an error
   u = (double *) malloc(n * sizeof(double));
   if (u == NULL) return -1;
           // For each column use a Householder transformation 
           //   to zero all entries below the subdiagonal.
   for (psubdiag = A + n, col = 0; col < (n - 2); psubdiag += (n+1), col++) {
          // Calculate the signed square root of the sum of squares of the
          // elements below the diagonal.
      for (pA = psubdiag, sss = 0.0, i = col + 1; i < n; pA += n, i++)
         sss += *pA * *pA;
      if (sss == 0.0) continue;
      sss = sqrt(sss);
      if ( *psubdiag >= 0.0 ) sss = -sss;
          // Calculate the Householder transformation Q = I - 2uu'/u'u.
      u[col + 1] = *psubdiag - sss;
      *psubdiag = sss; 
      for (pA = psubdiag + n, i = col + 2; i < n; pA += n, i++) {
         u[i] = *pA;
         *pA = 0.0;
      }
                    // Premultiply A by Q
      scale = -1.0 / (sss * u[col+1]);
      for (p_row = psubdiag - col, i = col + 1; i < n; i++) {
         pA = A + n * (col + 1) + i;
         for (innerproduct = 0.0, k = col + 1; k < n; pA += n, k++) 
            innerproduct += u[k] * *pA;
         innerproduct *= scale;
         for (pA = p_row + i, k = col + 1; k < n; pA += n, k++) 
            *pA -= u[k] * innerproduct;
      }
                    // Postmultiply QA by Q
      for (p_row = A, i = 0; i < n; p_row += n, i++) {
         for (innerproduct = 0.0, k = col + 1; k < n; k++) 
            innerproduct += u[k] * *(p_row + k);
         innerproduct *= scale;
         for (k = col + 1; k < n; k++) 
            *(p_row + k) -= u[k] * innerproduct;
      }
                    // Postmultiply U by (I - 2uu')
      for (i = 0, pU = U; i < n; pU += n, i++) {
         for (innerproduct = 0.0, k = col + 1; k < n; k++) 
            innerproduct += u[k] * *(pU + k);
         innerproduct *= scale;
         for (k = col + 1; k < n; k++) 
            *(pU + k) -= u[k] * innerproduct;
      }
   }
   free(u);
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  int QR_Hessenberg_Matrix( double *H, double *S, double eigen_real[],      //
//                     double eigen_imag[], int n, int Max_Iter)   //
//                                                                            //
//  Description:                                                              //
//     This program calculates the eigenvalues and eigenvectors of a matrix   //
//     in Hessenberg form. This routine is adapted from the routine 'hql2'    //
//     appearing in 'Handbook for Automatic Computation, vol 2: Linear        //
//     Algebra' published by Springer Verlag (1971) edited by Wilkinson and   //
//     Reinsch, Contribution II/15 Eigenvectors of Real and Complex Matrices  //
//     by LR and QR Triangularizations by Peters and Wilkinson.               //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the real n x n matrix H in upper//
//            Hessenberg form.                                                //
//     double *S                                                              //
//            If H is the primary data matrix, the matrix S should be set     //
//            to the identity n x n identity matrix on input.  If H is        //
//            derived from an n x n matrix A, then S should be the            //
//            transformation matrix such that AS = SH.  On output, the i-th   //
//            column of S corresponds to the i-th eigenvalue if that eigen-   //
//            value is real and the i-th and (i+1)-st columns of S correspond //
//            to the i-th eigenvector with the real part in the i-th column   //
//            and positive imaginary part in the (i+1)-st column if that      //
//            eigenvalue is complex with positive imaginary part.  The        //
//            eigenvector corresponding to the eigenvalue with negative       //
//            imaginary part is the complex conjugate of the eigenvector      //
//            corresponding to the complex conjugate of the eigenvalue.       //
//            If on input, S was the identity matrix, then the columns are    //
//            the eigenvectors of H as described, if S was a transformation   //
//            matrix so that AS = SH, then the columns of S are the           //
//            eigenvectors of A as described.                                 //
//     double eigen_real[]                                                    //
//            Upon return, eigen_real[i] will contain the real part of the    //
//            i-th eigenvalue.                                                //
//     double eigen_imag[]                                                    //
//            Upon return, eigen_ima[i] will contain the imaginary part of    //
//            the i-th eigenvalue.                                            //
//     int    n                                                               //
//            The number of rows or columns of the upper Hessenberg matrix A. //
//     int    Max_Iter                                             //
//            The maximum number of iterations to try to find an eigenvalue   //
//            before quitting.                                                //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - Unable to find an eigenvalue within 'Max_Iter' //
//                  iterations.                                               //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define Max_Iter                                            //
//     double H[N][N], S[N][N], eigen_real[N], eigen_imag[N];                 //
//     int k;                                                                 //
//                                                                            //
//     (code to initialize H[N][N] and S[N][N])                               //
//     k = QR_Hessenberg_Matrix( (double*)H, (double*)S, eigen_real,          //
//                                      eigen_imag, N, Max_Iter);  //
//     if (k < 0) {printf("Failed"); exit(1);}                                //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int QR_Hessenberg_Matrix( double *H, double *S, double eigen_real[],
                          double eigen_imag[], int n, unsigned int Max_Iter)
{
   int i;
   int row;
   unsigned int iteration;
   int found_eigenvalue;
   double shift = 0.0;
   double* pH;
   for ( row = n - 1; row >= 0; row--) {
      found_eigenvalue = 0;
      for (iteration = 1; iteration <= Max_Iter; iteration++) {
                      // Search for small subdiagonal element
         for (i = row, pH = H + row * n; i > 0; i--, pH -= n)
            if (fabs(*(pH + i - 1 )) <= DBL_EPSILON *
                       ( fabs(*(pH - n + i - 1)) + fabs(*(pH + i)) ) ) break;
                // If the subdiagonal element on row "row" is small, then
                // that row element is an eigenvalue.  If the subdiagonal 
                // element on row "row-1" is small, then the eigenvalues
                // of the 2x2 diagonal block consisting rows "row-1" and
                // "row" are eigenvalues.  Otherwise perform a double QR
                // iteration.
         switch(row - i) {
            case 0: // One real eigenvalue
               One_Real_Eigenvalue(pH, eigen_real, eigen_imag, i, shift);
               found_eigenvalue = 1;
               break;
            case 1: // Either two real eigenvalues or a complex pair
               row--;
               Two_Eigenvalues(H, S, eigen_real, eigen_imag, n, row, shift);
               found_eigenvalue = 1;
               break;    
            default:
               Double_QR_Iteration(H, S, i, row, n, &shift, iteration);
         }
         if (found_eigenvalue) break;
      }
      if (iteration > Max_Iter) return -1;
   }
   BackSubstitution(H, eigen_real, eigen_imag, n);
   Calculate_Eigenvectors(H, S, eigen_real, eigen_imag, n);
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  static void One_Real_Eigenvalue( double Hrow[], double eigen_real[],      //
//                                double eigen_imag[], int row, double shift) //
//                                                                            //
//  Arguments:                                                                //
//     double Hrow[]                                                          //
//            Pointer to the row "row" of the matrix in Hessenberg form.      //
//     double eigen_real[]                                                    //
//            Array of the real parts of the eigenvalues.                     //
//     double eigen_imag[]                                                    //
//            Array of the imaginary parts of the eigenvalues.                //
//     int    row                                                             //
//            The row to which the pointer Hrow[] points of the matrix H.     //
//     double shift                                                           //
//            The cumulative exceptional shift of the diagonal elements of    //
//            the matrix H.                                                   //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void One_Real_Eigenvalue(double Hrow[], double eigen_real[],
                                   double eigen_imag[], int row, double shift)
{
   Hrow[row] += shift;      
   eigen_real[row] = Hrow[row];
   eigen_imag[row] = 0.0;
}
////////////////////////////////////////////////////////////////////////////////
//  static void Two_Eigenvalues( double *H, double *S, double eigen_real[],   //
//                         double eigen_imag[], int n, int row, double shift) //
//                                                                            //
//  Description:                                                              //
//     Given the 2x2 matrix A = (a[i][j]), the characteristic equation is:    //
//     x^2 - Tr(A) x + Det(A) = 0, where Tr(A) = a[0][0] + a[1][1] and        //
//     Det(A) = a[0][0] * a[1][1] - a[0][1] * a[1][0].                        //
//     The solution for the eigenvalues x are:                                //
//         x1 = (Tr(A) + sqrt( (Tr(A))^2 + 4 * a[0][1] * a[1][0] ) / 2 and    //
//         x2 = (Tr(A) - sqrt( (Tr(A))^2 + 4 * a[0][1] * a[1][0] ) / 2.       //
//     Let p = (a[0][0] - a[1][1]) / 2 and q = p^2 - a[0][1] * a[1][0], then  //
//         x1 = a[1][1] + p [+|-] sqrt(q) and x2 = a[0][0] + a[1][1] - x1.    //
//     Choose the sign [+|-] to be the sign of p.                             //
//     If q > 0.0, then both roots are real and the transformation            //
//                 | cos sin |    | a[0][0]  a[0][1] |    | cos -sin |        //
//                 |-sin cos |    | a[1][0]  a[1][1] |    | sin  cos |        //
//      where sin = a[1][0] / r, cos = ( p + sqrt(q) ) / r, where r > 0 is    //
//      determined sin^2 + cos^2 = 1 transforms the matrix A to an upper      //
//      triangular matrix with x1 the upper diagonal element and x2 the lower.//
//      If q < 0.0, then both roots form a complex conjugate pair.            //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double *S                                                              //
//            Pointer to the first element of the transformation matrix.      //
//     double eigen_real[]                                                    //
//            Array of the real parts of the eigenvalues.                     //
//     double eigen_imag[]                                                    //
//            Array of the imaginary parts of the eigenvalues.                //
//     int    n                                                               //
//            The dimensions of the matrix H and S.                           //
//     int    row                                                             //
//            The upper most row of the block diagonal 2 x 2 submatrix of H.  //
//     double shift                                                           //
//            The cumulative exceptional shift of the diagonal elements of    //
//            the matrix H.                                                   //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Two_Eigenvalues(double *H, double* S, double eigen_real[],
                            double eigen_imag[], int n, int row, double shift)
{
   double p, q, x, discriminant, r;
   double cos, sin;
   double *Hrow = H + n * row;
   double *Hnextrow = Hrow + n;
   int nextrow = row + 1;
   p = 0.5 * (Hrow[row] - Hnextrow[nextrow]);
   x = Hrow[nextrow] * Hnextrow[row];
   discriminant = p * p + x;
   Hrow[row] += shift;
   Hnextrow[nextrow] += shift;
   if (discriminant > 0.0) {                 // pair of real roots
      q = sqrt(discriminant);
      if (p < 0.0) q = p - q; else q += p;
      eigen_real[row] = Hnextrow[nextrow] + q;
      eigen_real[nextrow] = Hnextrow[nextrow] - x / q;
      eigen_imag[row] = 0.0;
      eigen_imag[nextrow] = 0.0;
      r = sqrt(Hnextrow[row]*Hnextrow[row] + q * q);
      sin = Hnextrow[row] / r;
      cos = q / r;
      Update_Row(Hrow, cos, sin, n, row);
      Update_Column(H, cos, sin, n, row);
      Update_Transformation(S, cos, sin, n, row);
   }
   else {                             // pair of complex roots
      eigen_real[nextrow] = eigen_real[row] = Hnextrow[nextrow] + p;
      eigen_imag[row] = sqrt(fabs(discriminant));
      eigen_imag[nextrow] = -eigen_imag[row];
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Update_Row(double *Hrow, double cos, double sin, int n,       //
//                                                                   int row) //
//                                                                            //
//  Description:                                                              //
//     Update rows 'row' and 'row + 1' using the rotation matrix:             //
//                                | cos sin |                                 //
//                                |-sin cos |.                                //
//     I.e. multiply the matrix H on the left by the identity matrix with     //
//     the 2x2 diagonal block starting at row 'row' replaced by the above     //
//     2x2 rotation matrix.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double Hrow[]                                                          //
//            Pointer to the row "row" of the matrix in Hessenberg form.      //
//     double cos                                                             //
//            Cosine of the rotation angle.                                   //
//     double sin                                                             //
//            Sine of the rotation angle.                                     //
//     int    n                                                               //
//            The dimension of the matrix H.                                  //
//     int    row                                                             //
//            The row to which the pointer Hrow[] points of the matrix H      //
//            in Hessenberg form.                                             //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Update_Row(double *Hrow, double cos, double sin, int n, int row)
{
   double x;
   double *Hnextrow = Hrow + n;
   int i;
   for (i = row; i < n; i++) {
      x = Hrow[i];
      Hrow[i] = cos * x + sin * Hnextrow[i];
      Hnextrow[i] = cos * Hnextrow[i] - sin * x;
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Update_Column(double* H, double cos, double sin, int n,       //
//                                                                   int col) //
//                                                                            //
//  Description:                                                              //
//     Update columns 'col' and 'col + 1' using the rotation matrix:          //
//                               | cos -sin |                                 //
//                               | sin  cos |.                                //
//     I.e. multiply the matrix H on the right by the identity matrix with    //
//     the 2x2 diagonal block starting at row 'col' replaced by the above     //
//     2x2 rotation matrix.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the matrix in Hessenberg form.                       //
//     double cos                                                             //
//            Cosine of the rotation angle.                                   //
//     double sin                                                             //
//            Sine of the rotation angle.                                     //
//     int    n                                                               //
//            The dimension of the matrix H.                                  //
//     int    col                                                             //
//            The left-most column of the matrix H to update.                 //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Update_Column(double* H, double cos, double sin, int n, int col)
{
   double x;
   int i;
   int next_col = col + 1;
   for (i = 0; i <= next_col; i++, H += n) {
      x = H[col];
      H[col] = cos * x + sin * H[next_col];
      H[next_col] = cos * H[next_col] - sin * x;
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Update_Transformation(double *S, double cos, double sin,      //
//                                                             int n, int k)  //
//                                                                            //
//  Description:                                                              //
//     Update columns 'k' and 'k + 1' using the rotation matrix:              //
//                               | cos -sin |                                 //
//                               | sin  cos |.                                //
//     I.e. multiply the matrix S on the right by the identity matrix with    //
//     the 2x2 diagonal block starting at row 'k' replaced by the above       //
//     2x2 rotation matrix.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *S                                                              //
//            Pointer to the row "row" of the matrix in Hessenberg form.      //
//     double cos                                                             //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double sin                                                             //
//            Pointer to the first element of the transformation matrix.      //
//     int    n                                                               //
//            The dimensions of the matrix H and S.                           //
//     int    k                                                               //
//            The row to which the pointer Hrow[] points of the matrix H.     //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Update_Transformation(double *S, double cos, double sin,
                                                                 int n, int k)
{
   double x;
   int i;
   int k1 = k + 1;
   for (i = 0; i < n; i++, S += n) {
      x = S[k];
      S[k] = cos * x + sin * S[k1];
      S[k1] = cos * S[k1] - sin * x;
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Double_QR_Iteration(double *H, double *S, int min_row,        //
//               int max_row, int n, double* shift, unsigned int interation)  //
//                                                                            //
//  Description:                                                              //
//     Calculate the trace and determinant of the 2x2 matrix:                 //
//                        | H[k-1][k-1]  H[k-1][k] |                          //
//                        |  H[k][k-1]    H[k][k]  |                          //
//     unless iteration = 0 (mod 10) in which case increment the shift and    //
//     decrement the first k elements of the matrix H, then fudge the trace   //
//     and determinant by  trace = 3/2( |H[k][k-1]| + |H[k-1][k-2]| and       //
//     det = 4/9 trace^2.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the matrix H in Hessenberg form.                     //
//     double *S                                                              //
//            Pointer to the transformation matrix S.                         //
//     int    min_row                                                         //
//            The top-most row in which the off-diagonal element of H is      //
//            negligible.  If no such row exists, then min_row = 0.           //
//     int    max_row                                                         //
//            The maximum row of the block 2 x 2 diagonal matrix used to      //
//            estimate the two eigenvalues for the two implicit shifts.       //
//     int    n                                                               //
//            The dimensions of the matrix H and S.                           //
//     double *shift                                                          //
//            The cumulative exceptional shift of the diagonal elements of    //
//            the matrix H.                                                   //
//     int    iteration                                                       //
//            Current iteration count.                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Double_QR_Iteration(double *H, double *S, int min_row, int max_row,
                                          int n, double* shift, unsigned int iteration) 
{
   int k;
   double trace, det;
   Product_and_Sum_of_Shifts(H, n, max_row, shift, &trace, &det, iteration);
   k = Two_Consecutive_Small_Subdiagonal(H, min_row, max_row, n, trace, det);
   Double_QR_Step(H, min_row, max_row, k, trace, det, S, n);
}
////////////////////////////////////////////////////////////////////////////////
//  static void Product_and_Sum_of_Shifts(double *H, int n, int max_row,      //
//       double* shift, double *trace, double *det, unsigned int interation)  //
//                                                                            //
//  Description:                                                              //
//     Calculate the trace and determinant of the 2x2 matrix:                 //
//                        | H[k-1][k-1]  H[k-1][k] |                          //
//                        |  H[k][k-1]    H[k][k]  |                          //
//     unless iteration = 0 (mod 10) in which case increment the shift and    //
//     decrement the first k elements of the matrix H, then fudge the trace   //
//     and determinant by  trace = 3/2( |H[k][k-1]| + |H[k-1][k-2]| and       //
//     det = 4/9 trace^2.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the matrix H in Hessenberg form.                     //
//     int    n                                                               //
//            The dimension of the matrix H.                                  //
//     int    max_row                                                         //
//            The maximum row of the block 2 x 2 diagonal matrix used to      //
//            estimate the two eigenvalues for the two implicit shifts.       //
//     double *shift                                                          //
//            The cumulative exceptional shift of the diagonal elements of    //
//            the matrix H.  Modified if an exceptional shift occurs.         //
//     double *trace                                                          //
//            Returns the trace of the 2 x 2 block diagonal matrix starting   //
//            at the row/column max_row-1.  For an exceptional shift, the     //
//            trace is set as described above.                                //
//     double *det                                                            //
//            Returns the determinant of the 2 x 2 block diagonal matrix      //
//            starting at the row/column max_row-1.  For an exceptional shift,//
//            the determinant is set as described above.                      //
//     int    iteration                                                       //
//            Current iteration count.                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Product_and_Sum_of_Shifts(double *H, int n, int max_row,
                     double* shift, double *trace, double *det, unsigned int iteration) 
{
   double *pH = H + max_row * n;
   double *p_aux;
   int i;
   int min_col = max_row - 1;
   if ( ( iteration % 10 ) == 0 ) {
      *shift += pH[max_row];
      for (i = 0, p_aux = H; i <= max_row; p_aux += n, i++)
         p_aux[i] -= pH[max_row];
      p_aux = pH - n;
      *trace = fabs(pH[min_col]) + fabs(p_aux[min_col - 1]);
      *det = *trace * *trace;
      *trace *= 1.5;
   }
   else {
      p_aux = pH - n;
      *trace = p_aux[min_col] + pH[max_row];
      *det = p_aux[min_col] * pH[max_row] - p_aux[max_row] * pH[min_col];
   }
};
////////////////////////////////////////////////////////////////////////////////
//  static int Two_Consecutive_Small_Subdiagonal(double* H, int min_row,      //
//                              int max_row, int n, double trace, double det) //
//                                                                            //
//  Description:                                                              //
//     To reduce the amount of computation in Francis' double QR step search  //
//     for two consecutive small subdiagonal elements from row nn to row m,   //
//     where m < nn.                                                          //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     int    min_row                                                         //
//            The row in which to end the search (search is from upwards).    //
//     int    max_row                                                         //
//            The row in which to begin the search.                           //
//     int    n                                                               //
//            The dimension of H.                                             //
//     double trace                                                           //
//            The trace of the lower 2 x 2 block diagonal matrix.             //
//     double det                                                             //
//            The determinant of the lower 2 x 2 block diagonal matrix.       //
//                                                                            //
//  Return Value:                                                             //
//     Row with negligible subdiagonal element or min_row if none found.      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int Two_Consecutive_Small_Subdiagonal(double* H, int min_row,
                                 int max_row, int n, double trace, double det)
{
   double x, y ,z, s;
   double* pH;
   int i, k;
   for (k = max_row - 2, pH = H + k * n; k >= min_row; pH -= n, k--) {
      x = (pH[k] * ( pH[k] - trace ) + det) / pH[n+k] + pH[k+1];
      y = pH[k] + pH[n+k+1] - trace;
      z = pH[n + n + k + 1];
      s = fabs(x) + fabs(y) + fabs(z);
      x /= s;
      y /= s;
      z /= s;
      if (k == min_row) break;
      if ( (fabs(pH[k-1]) * (fabs(y) + fabs(z)) ) <= 
          DBL_EPSILON * fabs(x) *
             (fabs(pH[k-1-n]) + fabs(pH[k]) + fabs(pH[n + k + 1])) ) break; 
   }
   for (i = k+2, pH = H + i * n; i <= max_row; pH += n, i++) pH[i-2] = 0.0;
   for (i = k+3, pH = H + i * n; i <= max_row; pH += n, i++) pH[i-3] = 0.0;
   return k;
};
////////////////////////////////////////////////////////////////////////////////
//  static void Double_QR_Step(double *H, int min_row, int max_row,           //
//                                            int min_col, double *S, int n)  //
//                                                                            //
//  Description:                                                              //
//     Perform Francis' double QR step from rows 'min_row' to 'max_row'       //
//     and columns 'min_col' to 'max_row'.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     int    min_row                                                         //
//            The row in which to begin.                                      //
//     int    max_row                                                         //
//            The row in which to end.                                        //
//     int    min_col                                                         //
//            The column in which to begin.                                   //
//     double trace                                                           //
//            The trace of the lower 2 x 2 block diagonal matrix.             //
//     double det                                                             //
//            The determinant of the lower 2 x 2 block diagonal matrix.       //
//     double *S                                                              //
//            Pointer to the first element of the transformation matrix.      //
//     int    n                                                               //
//            The dimensions of H and S.                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Double_QR_Step(double *H, int min_row, int max_row, int min_col,
                                   double trace, double det, double *S, int n)
{
   double s, x, y, z;
   double a, b, c;
   double *pH;
   double *tH;
   double *pS;
   int i,j,k;
   int last_test_row_col = max_row - 1;
   k = min_col;
   pH = H + min_col * n;
   a = (pH[k] * ( pH[k] - trace ) + det) / pH[n+k] + pH[k+1];
   b = pH[k] + pH[n+k+1] - trace;
   c = pH[n + n + k + 1];
   s = fabs(a) + fabs(b) + fabs(c);
   a /= s;
   b /= s;
   c /= s;
   for (; k <= last_test_row_col; k++, pH += n) {
      if ( k > min_col ) {
         c = (k == last_test_row_col) ? 0.0 : pH[n + n + k - 1];
         x = fabs(pH[k-1]) + fabs(pH[n + k - 1]) + fabs(c);
         if ( x == 0.0 ) continue;
         a = pH[k - 1] / x;
         b = pH[n + k - 1] / x;
         c /= x;
      }
      s = sqrt( a * a + b * b + c * c );
      if (a < 0.0) s = -s;
      if ( k > min_col ) pH[k-1] = -s * x;
      else if (min_row != min_col) pH[k-1] = -pH[k-1];
      a += s;
      x = a / s;
      y = b / s;
      z = c / s;
      b /= a;
      c /= a;
                          // Update rows k, k+1, k+2
      for (j = k; j < n; j++) {
         a = pH[j] + b * pH[n+j];
         if ( k != last_test_row_col ) {
            a += c * pH[n + n + j];
            pH[n + n + j] -= a * z;
         }
         pH[n + j] -= a * y;
         pH[j] -= a * x;
      }
                             // Update column k+1
      j = k + 3;
      if (j > max_row) j = max_row;
      for (i = 0, tH = H; i <= j; i++, tH += n) {
         a = x * tH[k] + y * tH[k+1];
         if ( k != last_test_row_col ) {
            a += z * tH[k+2];
            tH[k+2] -= a * c;
         }
         tH[k+1] -= a * b;
         tH[k] -= a;           
      }
                        // Update transformation matrix
      for (i = 0, pS = S; i < n; pS += n, i++) {
         a = x * pS[k] + y * pS[k+1];
         if ( k != last_test_row_col ) {
            a += z * pS[k+2];
            pS[k+2] -= a * c;
         }
         pS[k+1] -= a * b;
         pS[k] -= a;
      }
   }; 
} 
////////////////////////////////////////////////////////////////////////////////
//  static void BackSubstitution(double *H, double eigen_real[],              //
//                                               double eigen_imag[], int n)  //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double eigen_real[]                                                    //
//            The real part of an eigenvalue.                                 //
//     double eigen_imag[]                                                    //
//            The imaginary part of an eigenvalue.                            //
//     int    n                                                               //
//            The dimension of H, eigen_real, and eigen_imag.                 //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void BackSubstitution(double *H, double eigen_real[],
                                                   double eigen_imag[], int n)
{
   double zero_tolerance;
   double *pH;
   int i, j, row;
                        // Calculate the zero tolerance
   pH = H;
   zero_tolerance = fabs(pH[0]);
   for (pH += n, i = 1; i < n; pH += n, i++)
      for (j = i-1; j < n; j++) zero_tolerance += fabs(pH[j]);
   zero_tolerance *= DBL_EPSILON;
                           // Start Backsubstitution
   for (row = n-1; row >= 0; row--) {
      if (eigen_imag[row] == 0.0) 
         BackSubstitute_Real_Vector(H, eigen_real, eigen_imag, row,
                                                           zero_tolerance, n);
      else if ( eigen_imag[row] < 0.0 ) 
         BackSubstitute_Complex_Vector(H, eigen_real, eigen_imag, row,
                                                           zero_tolerance, n);
   } 
}
////////////////////////////////////////////////////////////////////////////////
//  static void BackSubstitute_Real_Vector(double *H, double eigen_real[],    //
//             double eigen_imag[], int row,  double zero_tolerance, int n)   //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double eigen_real[]                                                    //
//            The real part of an eigenvalue.                                 //
//     double eigen_imag[]                                                    //
//            The imaginary part of an eigenvalue.                            //
//     int    row                                                             //
//     double zero_tolerance                                                  //
//            Zero substitute. To avoid dividing by zero.                     //
//     int    n                                                               //
//            The dimension of H, eigen_real, and eigen_imag.                 //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void BackSubstitute_Real_Vector(double *H, double eigen_real[],
                   double eigen_imag[], int row,  double zero_tolerance, int n)
{
   double *pH;
   double *pV;
   double x;
   double u[4];
   double v[2];
   int i,j,k;
   k = row;
   pH = H + row * n;
   pH[row] = 1.0;
   for (i = row - 1, pH -= n; i >= 0; i--, pH -= n) {
      u[0] = pH[i] - eigen_real[row];
      v[0] = pH[row];
      pV = H + n * k;
      for (j = k; j < row; j++, pV += n) v[0] += pH[j] * pV[row];
      if ( eigen_imag[i] < 0.0 ) {
         u[3] = u[0];
         v[1] = v[0];
      } else {
         k = i;
         if (eigen_imag[i] == 0.0) {
            if (u[0] != 0.0) pH[row] = - v[0] / u[0];
            else pH[row] = - v[0] / zero_tolerance;
         } else {
            u[1] = pH[i+1];
            u[2] = pH[n+i];
            x = (eigen_real[i] - eigen_real[row]);
            x *= x;
            x += eigen_imag[i] * eigen_imag[i];
            pH[row] = (u[1] * v[1] - u[3] * v[0]) / x; 
            if ( fabs(u[1]) > fabs(u[3]) )
               pH[n+row] = -(v[0] + u[0] * pH[row]) / u[1];
            else 
               pH[n+row] = -(v[1] + u[2] * pH[row]) / u[3];
         }
      }
   }    
}
////////////////////////////////////////////////////////////////////////////////
//  static void BackSubstitute_Complex_Vector(double *H, double eigen_real[], //
//              double eigen_imag[], int row,  double zero_tolerance, int n)  //
//                                                                            //
//  Description:                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double eigen_real[]                                                    //
//            The real part of an eigenvalue.                                 //
//     double eigen_imag[]                                                    //
//            The imaginary part of an eigenvalue.                            //
//     int    row                                                             //
//     double zero_tolerance                                                  //
//            Zero substitute. To avoid dividing by zero.                     //
//     int    n                                                               //
//            The dimension of H, eigen_real, and eigen_imag.                 //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void BackSubstitute_Complex_Vector(double *H, double eigen_real[],
                   double eigen_imag[], int row,  double zero_tolerance, int n)
{
   double *pH;
   double *pV;
   double x,y;
   double u[4];
   double v[2];
   double w[2];
   int i,j,k;
   k = row - 1;
   pH = H + n * row;
   if ( fabs(pH[k]) > fabs(pH[row-n]) ) {
      pH[k-n] = - (pH[row] - eigen_real[row]) / pH[k];
      pH[row-n] = -eigen_imag[row] / pH[k];
   }
   else 
      Complex_Division(-pH[row-n], 0.0,
           pH[k-n]-eigen_real[row], eigen_imag[row], &pH[k-n], &pH[row-n]);
   pH[k] = 1.0;
   pH[row] = 0.0;
   for (i = row - 2, pH = H + n * i; i >= 0; pH -= n, i--) {
      u[0] = pH[i] - eigen_real[row];
      w[0] = pH[row];
      w[1] = 0.0;
      pV = H + k * n;
      for (j = k; j < row; j++, pV+=n) {
         w[0] += pH[j] * pV[row - 1];
         w[1] += pH[j] * pV[row];
      }
      if (eigen_imag[i] < 0.0) {
         u[3] = u[0];
         v[0] = w[0];
         v[1] = w[1];
      } else {
         k = i;
         if (eigen_imag[i] == 0.0) {
            Complex_Division(-w[0], -w[1], u[0], eigen_imag[row], &pH[row-1],
                                                                    &pH[row]);
         }
         else {
            u[1] = pH[i+1];
            u[2] = pH[n + i];
            x = eigen_real[i] - eigen_real[row];
            y = 2.0 * x * eigen_imag[row];
            x = x * x + eigen_imag[i] * eigen_imag[i] 
                                          - eigen_imag[row] * eigen_imag[row];
            if ( x == 0.0 && y == 0.0 ) 
              x = zero_tolerance * ( fabs(u[0]) + fabs(u[1]) + fabs(u[2])
                                       + fabs(u[3]) + fabs(eigen_imag[row]) );
            Complex_Division(u[1]*v[0] - u[3] * w[0] + w[1] * eigen_imag[row],
                           u[1] * v[1] - u[3] * w[1] - w[0] * eigen_imag[row],
                                                   x, y, &pH[row-1], &pH[row]);
            if ( fabs(u[1]) > (fabs(u[3]) + fabs(eigen_imag[row])) ) {
               pH[n+row-1] = -w[0] - u[0] * pH[row-1]
                                           + eigen_imag[row] * pH[row] / u[1];
               pH[n+row] = -w[1] - u[0] * pH[row]
                                          - eigen_imag[row] * pH[row-1] / u[1];
            }
            else {
               Complex_Division(-v[0] - u[2] * pH[row-1], -v[1] - u[2]*pH[row],
                              u[3], eigen_imag[row], &pH[n+row-1], &pH[n+row]);
            }
         } 
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Calculate_Eigenvectors(double *H, double *S,                  //
//                          double eigen_real[], double eigen_imag[], int n)  //
//                                                                            //
//  Description:                                                              //
//     Multiply by transformation matrix.                                     //
//                                                                            //
//  Arguments:                                                                //
//     double *H                                                              //
//            Pointer to the first element of the matrix in Hessenberg form.  //
//     double *S                                                              //
//            Pointer to the first element of the transformation matrix.      //
//     double eigen_real[]                                                    //
//            The real part of an eigenvalue.                                 //
//     double eigen_imag[]                                                    //
//            The imaginary part of an eigenvalue.                            //
//     int    n                                                               //
//            The dimension of H, S, eigen_real, and eigen_imag.              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Calculate_Eigenvectors(double *H, double *S, double eigen_real[],
                                                   double eigen_imag[], int n)
{
   double* pH;
   double* pS;
   double x,y;
   int i,j,k;
   for (k = n-1; k >= 0; k--) {
      if (eigen_imag[k] < 0.0) {
         for (i = 0, pS = S; i < n; pS += n, i++) {
            x = 0.0;
            y = 0.0;
            for (j = 0, pH = H; j <= k; pH += n, j++) {
               x += pS[j] * pH[k-1];
               y += pS[j] * pH[k];
            }
            pS[k-1] = x;
            pS[k] = y;
         }
      } else if (eigen_imag[k] == 0.0) { 
         for (i = 0, pS = S; i < n; i++, pS += n) {
            x = 0.0;
            for (j = 0, pH = H; j <= k; j++, pH += n)
               x += pS[j] * pH[k];
            pS[k] = x;
        }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
//  static void Complex_Division(double x, double y, double u, double v,      //
//                                                    double* a, double* b)   //
//                                                                            //
//  Description:                                                              //
//    a + i b = (x + iy) / (u + iv)                                           //
//            = (x * u + y * v) / r^2 + i (y * u - x * v) / r^2,              //
//    where r^2 = u^2 + v^2.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//            Real part of the numerator.                                     //
//     double y                                                               //
//            Imaginary part of the numerator.                                //
//     double u                                                               //
//            Real part of the denominator.                                   //
//     double v                                                               //
//            Imaginary part of the denominator.                              //
//     double *a                                                              //
//            Real part of the quotient.                                      //
//     double *b                                                              //
//            Imaginary part of the quotient.                                 //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Complex_Division(double x, double y, double u, double v,
                                                         double* a, double* b)
{
   double q = u*u + v*v;
   *a = (x * u + y * v) / q;
   *b = (y * u - x * v) / q;
}

int main ()
{
   int max_iter {1000};
   //double tolerance {0.00001};
   //double omega {1.25};
   //(void) omega;
   std::cout << std::setprecision(5);   
   
   int num_rows {0}, num_cols {0};
   matrix<double> A;
   ReadMatrixMarketMatrix(A,num_rows,num_cols);
   vector <double> B;
   ReadMatrixMarketVector(B, num_cols);      
   vector<double> X = InitializeVectorToZero(num_rows);
   vector<double> array_A = InitializeVectorToZero(num_rows*num_cols);
   //double *array_A;
   //array_A = new double[num_rows*num_cols];
   // Allocate memory for U and V on the heap for efficiency
   double* array_U = new double[num_rows * num_cols];
   double* array_V = new double[num_rows * num_cols];
   // Copy A matrix into U and V arrays
   for (int i = 0; i < num_rows; ++i) {
     for (int j = 0; j < num_cols; ++j) {
       array_A[i * num_rows + j] = A[i][j];
       array_U[i * num_rows + j] = A[i][j];
       array_V[i * num_rows + j] = A[i][j];
     }
   }
   // Allocate memory for D, dummy_array, and X on the heap
   vector<double> D(num_cols);
   vector<double> dummy_array(num_cols);
   
   ReadMatrixMarketVector(B, num_cols);  
   //   PrintVector(X);
   //   WriteMatrixMarketVector(X,num_rows,num_cols);
   //   WriteMatrixMarketMatrix(A,num_rows,num_cols);
   
   // Using time point and system_clock
   std::chrono::time_point<std::chrono::system_clock> start, end;
 
   start = std::chrono::system_clock::now();	 
   // Begin to place new code here.
   QR_Input(&array_A[0], &B[0], &X[0], num_rows, num_cols, max_iter);  
   
   // End the new code here. 
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);
 
   std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
   
   // Memory allocation released here.         
   delete[] array_U;
   delete[] array_V;  
   array_A.clear();
   array_A.shrink_to_fit();  	
   A.clear();
   A.shrink_to_fit();
   B.clear();
   B.shrink_to_fit();
   X.clear();
   X.shrink_to_fit();
   D.clear();
   D.shrink_to_fit();
   dummy_array.clear();
   dummy_array.shrink_to_fit();
   return 0;
}
