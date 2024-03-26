#include <string.h>                                 // required for memcpy()
#include <iostream>
#include <math.h>		// required for fabs() and for sqrt() required for frexp() and ldexp()
#include <float.h>                  // required for DBL_EPSILON required for DBL_MIN_EXP and DBL_MAX_EXP
#include <fstream>
#include <vector>
#include <iomanip>
#define SUCCESS 0
#define UNDERFLOWS 1
#define OVERFLOWS 2
using namespace std;

template<typename T>
using matrix = std::vector<std::vector<T>>;

std::vector<double> CreateZeroVector(int num_rows);
int Output_Inverse_Power_Results (double array_A[], int num_cols, double X0[],double tolerance, int max_tries);
int Inverse_Power_Solve(double array_A[], int num_cols, double eigenvalue, double X0[],double tolerance, int max_tries);
int Subtract_Scalar_from_Diagonal(double array_A[], double x_scalar, int num_rows, int num_cols);
void Copy_Vector(double *d, double *s, int num_cols);
int Crout_LU_Decomposition_with_Pivoting(double *A, int pivot[], int num_cols);
int Crout_LU_with_Pivoting_Solve(double *LU, double B[], int pivot[], double x[], int num_cols);
int Equilibrate_Matrix(double *A, int num_rows, int num_cols, double r[], double c[]);
int Equilibrate_Right_Hand_Side(double B[], double r[], int num_cols);
int Unequilibrate_Solution(double x[], double c[], int num_cols);
int Multiply_Matrix_with_Vector(double X[], double array_A[], int num_rows, int num_cols, double X0[]);
int Eigenvalue_Inverse_Power_Method(double array_A[], int num_cols, double eigenvalue, double X0[], double tolerance, int max_tries);
//                    Required Externally Defined Routines
// Place Externally Define Routines here.
double Inner_Product(double u[], double v[], int n);
double Vector_Max_Norm(double v[], int n);
double Vector_L2_Norm(double v[], int n);
void   Divide_Vector_by_Scalar(double v[], double x, int n);
// Place my created externnaly defined routines here.
void PrintVector2D(const vector<double> A[], int num_rows);
void PrintVector(const double B[], int num_rows);

int ReadMatrixMarketFile(matrix<double>& A, int& num_rows, int& num_cols);
int ReadVectorMarketFile(vector<double>& B, int& num_cols);
int Inverse_Power_Input();

std::vector<double> CreateZeroVector(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
}

int Output_Inverse_Power_Results (double array_A[], int num_cols, double X0[],double tolerance, int max_tries)
{
  int i,j;
  int num_rows = num_cols;
  std::cout << std::setprecision (5) << endl;
  printf ("The tolerance is ");
  printf ("%6.4lf \n", tolerance);
  printf ("The number of max_tries is %d \n", max_tries);
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where the matrix A = \n");
  for (i = 0; i < num_rows; i++)
  {
    for (j = 0; j < num_cols; j++)
	{
	  printf ("%6.4f   ", array_A[i * num_rows + j]);
	}
    printf ("\n");
  }
  printf ("\n");
  FILE *myFile3;
  myFile3 = fopen("C.dat","w+");
  if (myFile3 == NULL)
  {
    printf("Error writing to file.\n");
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, num_rows, num_cols);
  std::cout << "The normalized eigenvector is X0 = " << endl; 
  PrintVector(&X0[0], num_cols); 
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << num_rows << " " << num_rows;
    outputFile << endl;
    for (i = 0; i < num_cols; i++)
      outputFile <<  1 << " " << i+1 << " " << X0[i] << endl;
    // Close the file
    outputFile.close();
  } 
  else {
    std::cout << "Error writing to file." << std::endl;
  }
  return 0;
}
int Inverse_Power_Solve(double array_A[], int num_cols, double eigenvalue, double X0[],double tolerance, int max_tries)
{
  int err;
// Begin your code here.
  err = Eigenvalue_Inverse_Power_Method(&array_A[0], num_cols, eigenvalue, &X0[0], tolerance, max_tries);
//  err = Eigenvalue_Rayleigh_Method(&array_A[0], num_cols, eigenvalue, &X[0], &X0[0], tolerance, max_tries);
  if (err == -1)
    printf(" The iteration failed to converge\n");
  if (err == -2)
    printf(" The iteration vector vanished\n");
  if (err == -3)
  {
    printf(" The estimate of the dominant eigenvalue vanished\n");
  }
// End your code here.
  Output_Inverse_Power_Results (&array_A[0], num_cols, &X0[0], tolerance, max_tries);
  return 0;
}
int Subtract_Scalar_from_Diagonal(double array_A[], double x_scalar, int num_rows, int num_cols)
{
  int i,n;
  n = ((num_rows < num_cols) ? num_rows : num_cols);
  for (i=0;i<n;i++)
  {
    array_A[i*num_cols+i] = array_A[i*num_cols+i] - x_scalar;
  }
  return array_A[0];
}
////////////////////////////////////////////////////////////////////////////////
// File: copy_vector.c                                                        //
// Routine(s):                                                                //
//    Copy_Vector                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Copy_Vector(double *d, double *s, int n)                             //
//                                                                            //
//  Description:                                                              //
//     Copy the n dimensional vector s(source) to the n dimensional           //
//     vector d(destination).  The memory locations of the source and         //
//     destination vectors must not overlap, otherwise the results            //
//     are installation dependent.                                            //
//                                                                            //
//  Arguments:                                                                //
//      double *d  Pointer to the first element of the destination vector d.  //
//      double *s  Pointer to the first element of the source vector s.       //
//      int    n   The number of elements of the source / destination vectors.//
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N],  vd[N];                                                   //
//                                                                            //
//     (your code to initialize the vector v)                                 //
//                                                                            //
//     Copy_Vector(vd, v, N);                                                 //
//     printf(" Vector vd is \n");                                            //
////////////////////////////////////////////////////////////////////////////////
void Copy_Vector(double *d, double *s, int num_cols)
{
   memcpy(d, s, sizeof(double) * num_cols);
}
////////////////////////////////////////////////////////////////////////////////
// File: crout_pivot.c                                                        //
// Routines:                                                                  //
//    Crout_LU_Decomposition_with_Pivoting                                    //
//    Crout_LU_with_Pivoting_Solve                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n)   //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to decompose a row interchanged       //
//     version of the n x n matrix A into a lower triangular matrix L and a   //
//     unit upper triangular matrix U such that A = LU.                       //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Crout_LU_with_Pivoting_Solve.                         //
//     (see below).                                                           //
//                                                                            //
//     The Crout method with partial pivoting is: Determine the pivot row and //
//     interchange the current row with the pivot row, then assuming that     //
//     row k is the current row, k = 0, ..., n - 1 evaluate in order the      //
//     the following pair of expressions                                      //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to intialize the matrix A)                                  //
//                                                                            //
//     err = Crout_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);        //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Crout_LU_Decomposition_with_Pivoting(double *A, int pivot[], int num_cols)
{
   int i, j, k;
   double *p_k, *p_row, *p_col;
   double max;
//         For each row and column, k = 0, ..., n-1,
   for (k = 0, p_k = A; k < num_cols; p_k += num_cols, k++) {
//            find the pivot row
      pivot[k] = k;
      max = fabs( *(p_k + k) );
      for (j = k + 1, p_row = p_k + num_cols; j < num_cols; j++, p_row += num_cols) {
         if ( max < fabs(*(p_row + k)) ) {
            max = fabs(*(p_row + k));
            pivot[k] = j;
            p_col = p_row;
         }
      }
//     and if the pivot row differs from the current row, then
//     interchange the two rows.
      if (pivot[k] != k)
         for (j = 0; j < num_cols; j++) {
            max = *(p_k + j);
            *(p_k + j) = *(p_col + j);
            *(p_col + j) = max;
         }
//                and if the matrix is singular, return error
      if ( *(p_k + k) == 0.0 ) return -1;
//      otherwise find the upper triangular matrix elements for row k.
      for (j = k+1; j < num_cols; j++) {
         *(p_k + j) /= *(p_k + k);
      }
//            update remaining matrix
      for (i = k+1, p_row = p_k + num_cols; i < num_cols; p_row += num_cols, i++)
         for (j = k+1; j < num_cols; j++)
            *(p_row + j) -= *(p_row + k) * *(p_k + j);
   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_with_Pivoting_Solve(double *LU, double B[], int pivot[],     //
//                                                        double x[], int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to solve the linear equation Ax = B.  //
//     This routine is called after the matrix A has been decomposed into a   //
//     product of a lower triangular matrix L and a unit upper triangular     //
//     matrix U without pivoting.  The argument LU is a pointer to the matrix //
//     the superdiagonal part of which is U and the subdiagonal together with //
//     the diagonal part is L. (The diagonal part of U is 1 and is not        //
//     stored.)   The matrix A = LU.                                          //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *LU      Pointer to the first element of the matrix whose       //
//                     elements form the lower and upper triangular matrix    //
//                     factors of A.                                          //
//     double *B       Pointer to the column vector, (n x 1) matrix, B.       //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     double *x       Solution to the equation Ax = B.                       //
//     int     n       The number of rows or columns of the matrix LU.        //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Crout_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);        //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else {                                                                 //
//        err = Crout_LU_with_Pivoting_Solve(&A[0][0], B, pivot, x, n);       //
//        if (err < 0) printf(" Matrix A is singular\n");                     //
//        else printf(" The solution is \n");                                 //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Crout_LU_with_Pivoting_Solve(double *LU, double B[], int pivot[], double x[], int num_cols)
{
   int i, k;
   double *p_k;
   double dum;
//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.
   for (k = 0, p_k = LU; k < num_cols; p_k += num_cols, k++) {
      if (pivot[k] != k) {dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
      x[k] = B[k];
      for (i = 0; i < k; i++) x[k] -= x[i] * *(p_k + i);
      x[k] /= *(p_k + k);
   }
//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.
   for (k = num_cols-1, p_k = LU + num_cols*(num_cols-1); k >= 0; k--, p_k -= num_cols) {
      if (pivot[k] != k) {dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
      for (i = k + 1; i < num_cols; i++) x[k] -= x[i] * *(p_k + i);
      if (*(p_k + k) == 0.0) return -1;
   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
// File: equilibrate_matrix.c                                                 //
// Routines:                                                                  //
//    Equilibrate_Matrix                                                      //
//    Equilibrate_Right_Hand_Side                                             //
//    Unequilibrate_Solution                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The purpose for equilibration of a matrix is to transform the matrix   //
//     to another matrix with a lower condition number by scaling the         //
//     the rows and columns of the original matrix.                           //
//                                                                            //
//     Given a system of linear equations  Ax = B, where A is an nxn square   //
//     matrix and x and B are vectors (nx1 matrices), the matrix A is         //
//     equilibrated by finding real diagonal matrices R and C such that the   //
//     condition number of RAC is less than that of the matrix A.  There are  //
//     infinitely many possible algorithms for the row scaling factors and    //
//     column scaling factors, but the most common methods are:               //
//        (1)  Set r[i] = 1 / max(|a[i][j]| : j = 0,...,n-1), i=0,...,n-1.    //
//             Set c[j] = 1 / max(|r[i]a[i][j]| : i = 0,...,n-1), j=0,..,n-1. //
//        (2)  Set r[i] = 1 / Sum |a[i][j]|, where the sum is over j,         //
//                                                               i=0,...,n-1. //
//             Set c[j] = 1.0, j=0,...,n-1                                    //
//        (3)  A variation of method (1) in which                             //
//             r[i] = 1/2^k, k = max( m : |a[i][j]| = x 2^m, 1/2 < x <= 1,    //
//                                              j = 0,...,n-1), i=0,...,n-1.  //
//             c[j] = 1/2^k, k = max( m : |r[i]a[i][j]| = x 2^m, 1/2 < x <= 1,//
//                                              i = 0,...,n-1), j=0,...,n-1.  //
//             By choosing the scaling factors as powers of 2, there is no    //
//             loss in the number of significant bits of the scaled matrix    //
//             elements.  Note that if a[i][j] = 0, then m is not defined.    //
//             If the i-th row consists entirely of 0's, then r[i] is set to  //
//             1.0 and if the j-th column consists entirely of 0's, then c[j] //
//             is set to 1.0.                                                 //
//                                                                            //
//     The equilibrated system, RACy = RB, is then solved for y using a       //
//     linear equation solver.  The solution of the original system of linear //
//     equations, the unequilibrated system, is x = Cy.  The advantage of     //
//     method (2) is that C is the identity so x = y.                         //
//                                                                            //
//     In the routine below, Equilibrate_Matrix(), method (3) is used to      //
//     equilibrate the matrix A where the diagonal scaling factors are        //
//     returned in the array r[] and the column scaling factors are returned  //
//     in the array c[].  The routine, Equilibrate_Matrix(), returns 0 if     //
//     equilibration of the matrix A is successful, bit 0 is set if           //
//     equilibration results in an underflow, i.e. a loss of significance of  //
//     a matrix element of RAC.  The matrix A is transformed as A -> RAC,     //
//     with a possible underflow if the return bit 0 is set.                  //
//                                                                            //
//     The routine, Equilibrate_Right_Hand_Side(), the diagonal elements of   //
//     or the row scaling elements are used to transform the right-hand side  //
//     of the linear equation Ax=B. Equilibrate_Right_Hand_Side() returns a 0 //
//     if multiplication of R by the vector B is successful, bit 0 is set if  //
//     multiplication of R by B results in an underflow, or bit 1 is set if   //
//     multiplication of R by B results in an overflow.  If an overflow is    //
//     detected, then the multiplication process is terminated immediately.   //
//     The complex vector B is transformed as B -> RB.  In case of an         //
//     overflow, the vector B was modified but not completely transformed.    //
//                                                                            //
//     If the matrix A is successfully equilibrated and the vector B is       //
//     successfully transformed, then a linear equation solver can be called  //
//     to solve the equilibrated system RACy = RB for y.                      //
//                                                                            //
//     The routine, Unequilibrate_Solution(), is then used to transform the   //
//     solution of the equilibrated system to a solution of the original      //
//     unequilibrated system of linear equations.  Unequilibrate_Solution()   //
//     returns a 0 if multiplication of C by the solution vector of the       //
//     equilibrated system, y, is successful, bit 0 is set if multiplication  //
//     of C by y results in an underflow, or bit 1 is set if multiplication   //
//     of C by y results in an overflow.  If an overflow is detected, then the//
//     multiplication process is terminated immediately, otherwise the vector //
//     y is transformed as y -> Cy.                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// int Equilibrate_Matrix(double *A, int num_rows, int num_cols, double r[],        //
//                                                                double c[]) //
//                                                                            //
//  Description:                                                              //
//     This routine is called to equilibrate the elements of the possibly     //
//     ill-conditioned num_rows x num_cols matrix A.  Note that for solving a system//
//     of linear equations, the matrix A is a square matrix so that the number//
//     of rows and the number of columns are equal.                           //
//                                                                            //
//     The matrix A is equilibrated by finding diagonal matrices R and C such //
//     that the absolute value of each element of the equilibrated matrix RAC //
//     is bounded by 1.  By choosing the scaling factors as powers of 2, there//
//     is no loss in the number of significant bits of the scaled matrix      //
//     elements.                                                              //
//                                                                            //
//     In the routine below, Equilibrate_Matrix(), the diagonal elements of   //
//     the column scaling factors are returned in the array c[] while the     //
//     the diagonal elements of the row scaling factors are returned in the   //
//     in the array r[].  The routine, Equilibrate_Matrix(), returns a 0 if   //
//     equilibration of the matrix A is successful, bit 0 is set if           //
//     equilibration results in an underflow, i.e. a loss of significance of  //
//     a matrix element of A.  The matrix A is transformed as A -> RAC, with  //
//     a possible underflow if the return bit 0 is set.                       //
//                                                                            //
//  Arguments:                                                                //
//    double *A                                                               //
//       The address of the first element of the matrix A.                    //
//       Declared in the calling routine as double A[num_rows][num_cols].           //
//       On return the matrix A contains the equilibrated matrix.             //
//    int    num_rows                                                            //
//       The number of rows of the matrix A.                                  //
//    int    num_cols                                                            //
//       The number of columns of the matrix A.                               //
//    double r[]                                                              //
//       On return, the array containing the row scale factors.               //
//       Declared in the calling routine as r[num_rows].                         //
//    double c[]                                                              //
//       On return, the array containing the column scale factors.            //
//       Declared in the calling routine as c[num_cols].                         //
//                                                                            //
//    Return Values:                                                          //
//        0 Success                                                           //
//        1 Scaling results in an underflow, i.e. loss of significance.       //
//                                                                            //
//        The matrix A is transformed as A -> RAC, with a possible underflow  //
//        if the return bit 0 is set.                                         //
//                                                                            //
//        The array r[] is set with the diagonal of the diagonal matrix R,    //
//        and the array c[] is set with the diagonal of the diagonal matrix C.//
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[N][M], r[N], c[M];                                            //
//     int return_value;                                                      //
//                                                                            //
//     (your code to initialize the matrix A )                                //
//                                                                            //
//     return_value = Equilibrate_Matrix(&A[0][0], N, M, r, c);               //
//     ...                                                                    //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Equilibrate_Matrix(double *A, int num_rows, int num_cols, double r[], double c[])
{
   double *pA, *pC;
   double x;
   int min, max, exponent;
   int nonzero;
   int i,j;
   int return_code = SUCCESS;
                          // Find row scale factors.
   for (i = 0, pA = A; i < num_rows; i++, pA += num_cols) {
      max = DBL_MIN_EXP;
      min = DBL_MAX_EXP;
      nonzero = 0;
      r[i] = 1.0;
      for (j = 0, pC = pA; j < num_cols; j++, pC++) {
         x = frexp(*pC, &exponent);
         if ( x == 0.0 ) continue;
         if (exponent > max) max = exponent;
         if (exponent < min) min = exponent;
         nonzero = 1;
      }
      if (nonzero) {
         if ( min - max < DBL_MIN_EXP ) return_code |= UNDERFLOWS;
         r[i] = ldexp(1.0, -max);
         for (j = 0, pC = pA; j < num_cols; j++, pC++) *pC *= r[i];
      }
   }
                         // Find Column Scale Factors.
   for (i = 0, pA = A; i < num_cols; i++, pA++) {
      max = DBL_MIN_EXP;
      min = DBL_MAX_EXP;
      nonzero = 0;
      c[i] = 1.0;
      for (j = 0, pC = pA; j < num_rows; j++, pC += num_cols) {
         x = frexp(*pC, &exponent);
         if ( x == 0.0 ) continue;
         if (exponent > max) max = exponent;
         if (exponent < min) min = exponent;
         nonzero = 1;
      }
      if (nonzero) {
         if ( min - max < DBL_MIN_EXP ) return_code = UNDERFLOWS;
         c[i] = ldexp(1.0, -max);
         for (j = 0, pC = pA; j < num_rows; j++, pC += num_cols) *pC *= c[i];
      }
   }
   return return_code;
}
////////////////////////////////////////////////////////////////////////////////
//  int Equilibrate_Right_Hand_Side(double B[], double r[], int n)            //
//                                                                            //
//  Description:                                                              //
//     This routine multiplies each element of B[] by the corresponding       //
//     element of r[].  I.e. B[i] -> B[i]*r[i], for i = 0,...,n-1.            //
//                                                                            //
//  Arguments:                                                                //
//    double B[]                                                              //
//       The address of the first element of the vector B.                    //
//       Declared in the calling routine as double B[n].                      //
//       On return, B is transformed to RB.                                   //
//    double r[]                                                              //
//       The address of the first element of the vector r which was calculated//
//       by calling Equilibrate_Matrix().                                     //
//       Declared in the calling routine as double r[n].                      //
//    int    n                                                                //
//       The number of elements in the vectors B and r.                       //
//                                                                            //
//    Return Values:                                                          //
//        0 Success                                                           //
//        1 Scaling results in an underflow, i.e. loss of significance.       //
//        2 Scaling results in an overflow, process halted.                   //
//        3 Scaling results in both an underflow and an overflow, process     //
//          halted.                                                           //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], r[N], c[N];                                      //
//                                                                            //
//     (your code to initialize the vector B and matrix A)                    //
//                                                                            //
//     return_code = Equilibrate_Matrix(&A[0][0], N, N, r, c);                //
//     printf("return codes for equilibrate matrix\n");                       //
//     if (return_code == 0) printf("success\n");                             //
//     if (return_code &= 1) printf("warning loss of significance\n");        //
//     return_code = Equilibrate_Right_Hand_Side(B, r, N);                    //
//     printf("return codes for equilibrate right hand side\n");              //
//     if (return_code == 0) printf("success\n");                             //
//     if (return_code &= 1) printf("warning loss of significance\n");        //
//     if (return_code &= 2) printf("fatal error - overflow \n");             //
//     ...                                                                    //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Equilibrate_Right_Hand_Side(double B[], double r[], int num_cols)
{
   int i;
   int ret_code = SUCCESS;
   int B_exponent, r_exponent;
   double rx;
   (void) rx;
   double Bx;
   (void) Bx;
   for (i = 0; i < num_cols; i++) {
      Bx = frexp(B[i], &B_exponent);
      rx = frexp(r[i], &r_exponent);
      if (r_exponent + B_exponent < DBL_MIN_EXP) ret_code |= UNDERFLOWS;
      if (r_exponent + B_exponent > DBL_MAX_EXP) return ret_code += OVERFLOWS;
      B[i] *= r[i];
   }
   return ret_code;
}
////////////////////////////////////////////////////////////////////////////////
//  int Unequilibrate_Solution(double x[], double c[], int n)                 //
//                                                                            //
//  Description:                                                              //
//     This routine multiplies each element of x[] by the corresponding       //
//     element of c[].  I.e. x[i] -> x[i]*c[i], for i = 0,...,n-1.            //
//                                                                            //
//  Arguments:                                                                //
//    double *x                                                               //
//       The address of the first element of the vector x.                    //
//       Declared in the calling routine as double complex x[n].              //
//       On input, x is the solution of the equilibrated system of linear     //
//       equations, and on output, x is the solution of the unequilibrated    //
//       original system of linear equations provided the return code is 0.   //
//    double c[]                                                              //
//       The address of the first element of the vector c which was calculated//
//       by calling Equilibrate_Matrix().                                     //
//       Declared in the calling routine as double c[n].                      //
//    int    n                                                                //
//       The number of elements in the vectors x and c.                       //
//                                                                            //
//    Return Values:                                                          //
//        0 Success                                                           //
//        1 Scaling results in an underflow, i.e. loss of significance.       //
//        2 Backscaling results in an overflow, process halted.               //
//        3 Backscaling results both an underflow and an overflow, process    //
//          halted.                                                           //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double x[N];                                                           //
//     double A[N][N], B[N], r[N], c[N];                                      //
//                                                                            //
//     (your code to initialize the vector B, and matrix A)                   //
//                                                                            //
//     return_code = Equilibrate_Matrix(&A[0][0], N, N, r, c);                //
//     printf("return codes for equilibrate matrix\n");                       //
//     if (return_code == 0) printf("success\n");                             //
//     if (return_code &= 1) printf("warning loss of significance\n");        //
//     return_code = Equilibrate_Right_Hand_Side(B, r, N);                    //
//     printf("return codes for equilibrate right hand side\n");              //
//     if (return_code == 0) printf("success\n");                             //
//     if (return_code &= 1) printf("warning loss of significance\n");        //
//     if (return_code &= 2) printf("fatal error - overflow \n");             //
//     if (return_code &= 2) exit(0);                                         //
//     (Call linear equation solver, return solution in x.)                   //
//     (If return from equation solver fails then exit)                       //
//     return_code = Unequilibrate_Solution(x, c, N);                         //
//     printf("return codes for unequilibrate solution\n");                   //
//     if (return_code == 0) printf("success\n");                             //
//     if (return_code &= 1) printf("warning loss of significance\n");        //
//     if (return_code &= 2) printf("fatal error - overflow \n");             //
//     ...                                                                    //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Unequilibrate_Solution(double x[], double c[], int num_cols)
{
   int i;
   int ret_code = SUCCESS;
   int c_exponent, x_exponent;
   double xx;
   (void) xx;
   double cx;
   (void) cx;
   for (i = 0; i < num_cols; i++) {
      xx = frexp(x[i], &x_exponent);
      cx = frexp(c[i], &c_exponent);
      if (c_exponent + x_exponent < DBL_MIN_EXP) ret_code = UNDERFLOWS;
      if (c_exponent + x_exponent > DBL_MAX_EXP) return ret_code += OVERFLOWS;
      x[i] *= c[i];
   }
   return ret_code;
};
////////////////////////////////////////////////////////////////////////////////
// File: multiply_matrix_by_vector.c                                          //
// Routine(s):                                                                //
//    Multiply_Matrix_by_Vector                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Matrix_by_Vector(double u[], double *A, int num_rows,          //
//                                                   int num_cols, double v[])   //
//                                                                            //
//  Description:                                                              //
//     Post multiply the num_rows x num_cols matrix A by the column vector v        //
//     to form the column vector u, i.e. u = A v.                             //
//     The matrix A should be declared as "double A[num_rows][num_cols]" in the     //
//     calling routine.  The vector v declared as "double v[num_cols]" and       //
//     the vector u declared as "double u[num_rows]" in the calling routine.     //
//                                                                            //
//  Arguments:                                                                //
//     double *u    Pointer to the first element of the vector u.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    num_rows The number of rows of the matrix A and the number of      //
//                  components of the column vector u.                        //
//     int    num_cols The number of columns of the matrices A and the           //
//                  number of components of the column vector v.              //
//     double *v    Pointer to the first element of the vector v.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  u[M], v[N];                                           //
//                                                                            //
//     (your code to initialize the matrix A and column vector v)             //
//                                                                            //
//     Multiply_Matrix_by_Vector(u, &A[0][0], M, N, v);                       //
//     printf("The vector u is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Multiply_Matrix_by_Vector(double u[], double *A, int num_rows, int num_cols, double v[])
{
   int i,j;
   for (i = 0; i < num_rows; A += num_cols, i++)
      for (u[i] = 0.0, j = 0; j < num_cols; j++) u[i] += A[j] * v[j];
}
////////////////////////////////////////////////////////////////////////////////
//  int Eigenvalue_Inverse_Power_Method(double *A, int n, double *eigenvalue, //
//                              double x0[], double tolerance, int max_tries) //
//                                                                            //
//  Description:                                                              //
//     The inverse power method applies the power method to the inverse of    //
//     the matrix (A - zI).  Since the power method yields the dominant       //
//     eigenvalue, the inverse power method yields the eigenvalue nearest     //
//     z.  If m is the dominant eigenvalue of the inverse of (A - zI), then   //
//     z + 1/m is the eigenvalue of A nearest z.                              //
//     The inverse power method can therefore be used to estimate the eigen-  //
//     values of a matrix A other than the dominant eigenvalue.               //
//                                                                            //
//     The inverse power method is:  (1) Start with the initial estimate      //
//     of the eigenvalue, z, and take the LU decomposition of the             //
//     matrix A - zI, (2) Start with an initial non-zero vector x0, (3)       //
//     Normalize the current vector x[k], (4) Solve (A-zI)x[k+1] = x[k], for  //
//     x[k+1], (5) let m = x[k+1]'x[k], (5) if the relative difference of the //
//     value of m and the old value of m is less than a preassigned tolerance,//
//     then halt the procedure; otherwise go to step (3).  The estimate of the//
//     eigenvalue of A nearest z is z + 1/m.                                  //
//                                                                            //
//     Note! If the matrix A does not admit n linearly independent eigen-     //
//     vectors, a dominant eigenvalue may still be found.                     //
//                                                                            //
//     Note! If there are more than one dominant eigenvalue of the inverse of //
//     (A - zI) the procedure may not converge.  If the dominant eigenvalue   //
//     is complex, then the complex conjugate is also a dominant eigenvalue.  //
//     And some refinements to this procedure need to be made.                //
//                                                                            //
//     Note! If the initial vector x0 belongs to the orthogonal complement    //
//     of the eigenvector corresponding to the dominant eigenvalue, then      //
//     theoretically all subsequent iterates will also.  But round-off errors //
//     will contaminate this theoretical result and introduce non-zero        //
//     components of x[k] along the direction of the eigenvector corresponding//
//     to the dominant eigenvalue.  Convergence then behaves as before.       //
//     If, however, the initial vector x0 is an eigenvector corresponding to  //
//     eigenvalue of the inverse of (A-zI), then the procedure will terminate.//
//     In particular if the inverse power method is used to estimate several  //
//     eigenvalues, then do not use the values of x0 returned from a call to  //
//     Eigenvalues_Inverse_Power_Method as the input values to the next call. //
//                                                                            //
//  Arguments:                                                                //
//     double *A                                                              //
//        The pointer to the first element of the matrix A[n][n].  The matrix //
//        A is destroyed during the computations.                             //
//     int     n                                                              //
//        The number of rows and/or columns of the matrix A.                  //
//     double *eigenvalue                                                     //
//        On input, *eigenvalue contains the approximate eigenvalue.  On      //
//        output, if the call returns a positive number indicating success,   //
//        then the estimate of the eigenvalue nearest the input value         //
//        replaces that in the address pointed to by *eigenvalue.             //
//     double x0[]                                                            //
//        On input the initial vector.  It should be declared in the calling  //
//        routine as "double x0[n]", where n is as above.  It should be       //
//        initialized to a non-zero vector.                                   //
//     double tolerance                                                       //
//        The user specified tolerance controls the halting of the procedure. //
//        When two successive estimates of the dominant eigenvalue have a     //
//        relative absolute difference less than tolerance the iteration      //
//        halts. If the user specified tolerance is less than the machine     //
//        epsilon, the machine epsilon is used instead of the user specified  //
//        tolerance.                                                          //
//     int max_tries                                                          //
//        A control parameter which governs the maximum number of iterations  //
//        to try before giving up.                                            //
//                                                                            //
//  Return Values:                                                            //
//     n >=0  Success - Value returned is the number of iterations performed. //
//      -1    Failure - Failed to converge within max_tries iterations.       //
//      -2    Failure - Initial vector or subsequent vector is 0.             //
//      -3    Failure - Estimate of eigenvalue blows up.                      //
//      -4    Failure - Unable to allocate enough dynamic memory.             //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], x[N], x0[N];                                           //
//     double eigenvalue, tolerance;                                          //
//     int pivot[N], c[N], r[N];                                              //
//     int max_tries;                                                         //
//                                                                            //
//     (your code to initialize the matrix A, the vector x0, tolerance,       //
//      initial estimate of eigenvalue and max_tries. )                       //
//     err = Eigenvalue_Inverse_Power_Method((double*)A, N, &eigenvalue, x0,  //
//                                                     tolerance, max_tries); //
//     if (err == -1)                                                         //
//        printf(" The iteration failed to converge\n");                      //
//     else if (err == -2)                                                    //
//        printf(" The iteration vector vanished\n");                         //
//     else if (err == -3)                                                    //
//        printf(" The estimate of the dominant eigenvalue blows up\n");      //
//     else if (err == -4)                                                    //
//        printf(" Ran out of memory.\n");                                    //
//     else printf(" The eigenvalue is \n");                                  //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Eigenvalue_Inverse_Power_Method(double array_A[], int num_cols, double eigenvalue, double X0[], double tolerance, int max_tries)
{
   int    i, j, tries, error = 0;
   double old_estimate;
   double norm;
   double lambda;
   double *px0;
   double *pdum, *pc, *pr;
   double *x;
   int* pivot;
   double A[num_cols][num_cols];
   for (i = 0; i < num_cols; i++)
   {
     for (j = 0; j < num_cols; j++)
	 {
	   A[i][j] = array_A[i * num_cols + j];
	 }
   }
        // Check that the user specified tolerance is greater than the
        // machine epsilon.  If not, set it to the machine epsilon.
   if (tolerance < DBL_EPSILON) tolerance = DBL_EPSILON;
                // Check that the initial vector is non-zero.
   norm = Vector_L2_Norm(X0, num_cols);
   if ( norm == 0.0 )  return -2;
               // Dynamically allocate working storage memory.
   pivot = (int*) malloc( sizeof(int) * num_cols);
   x   = (double*) malloc( sizeof(double) * num_cols);
   px0 = (double*) malloc( sizeof(double) * num_cols);
   pc  = (double*) malloc( sizeof(double) * num_cols);
   pr  = (double*) malloc( sizeof(double) * num_cols);
   if (pr == NULL) { error = -4; goto ret; }
                       // Normalize the initial vector.
   Divide_Vector_by_Scalar(X0, norm, num_cols);
          // Subtract the approximate eigenvalue from the diagonal of A
          // and calculate the LU decomposition of A - *eigenvalue I.
          // Return immediately if *eigenvalue is already an eigenvalue.
   Subtract_Scalar_from_Diagonal(&array_A[0], eigenvalue, num_cols, num_cols);
   Equilibrate_Matrix(*A, num_cols, num_cols, pr, pc);
   if ( Crout_LU_Decomposition_with_Pivoting(*A, pivot, num_cols) < 0 ) goto ret;
           // Calculate initial estimate of the dominant eigenvalue.
   Copy_Vector(px0, X0, num_cols);
   Equilibrate_Right_Hand_Side(px0, pr, num_cols);
   Crout_LU_with_Pivoting_Solve(*A, px0, pivot, x, num_cols);
   Unequilibrate_Solution(x, pc, num_cols);
   lambda = Inner_Product(x, X0, num_cols);
      // Iterate for a maximum of max_tries times using the inverse power
      // method.  Exit if an error occurs or if the relative difference of
      // two successive estimates is less than the specified tolerance.
   for (tries = 1; tries < max_tries; tries++) {
      old_estimate = lambda;
      pdum = x;
      x = X0;
      X0 = pdum;
      norm = Vector_L2_Norm(X0, num_cols);
      if ( norm == 0.0 )  {error = -2; goto ret; }
      Divide_Vector_by_Scalar(X0, norm, num_cols);
      Copy_Vector(px0, X0, num_cols);
      Equilibrate_Right_Hand_Side(px0, pr, num_cols);
      Crout_LU_with_Pivoting_Solve(*A, px0, pivot, x, num_cols);
      Unequilibrate_Solution(x, pc, num_cols);
      lambda = Inner_Product(x, X0, num_cols);
      if (lambda == 0.0) {error = -3; goto ret; }
      if (fabs((lambda - old_estimate) / lambda) < tolerance)
      {
        double revised_lambda = old_estimate;
        double dominant_eigenvalue = 1.0/revised_lambda;
        printf("The dominant eigenvalue is %6.4lf .\n", dominant_eigenvalue);
        break;
      }
   }
   eigenvalue += 1.0 / lambda;
   if (tries >= max_tries) error = -1;
ret:
   free(pr);
   free(pc);
   free(px0);
   free(x);
   free(pivot);
   if (error != 0) return error;
   return tries;
}
////////////////////////////////////////////////////////////////////////////////
// File: inner_product.c                                                      //
// Routines:                                                                  //
//    Inner_Product                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Inner_Product(double u[], double v[], int n)                       //
//                                                                            //
//  Description:                                                              //
//     An inner product on a real vector space V is a positive definite       //
//     symmetric bilinear form < , > : V x V -> R.                            //
//     Given a Real vector space V with a basis e[i] for which the inner      //
//     product of the basis vectors <e[i],e[j]> = delta(i,j), delta(i,j) being//
//     the Kronecker delta function, the inner product of two vectors in V    //
//     is the sum of the component-wise products.  If dim V = 3,              //
//     u = u[0] e[0] + ... + u[n-1] e[n-1] and                                //
//     v = v[0] e[0] + ... + v[n-1] e[n-1] are vectors in V, then             //
//     <u,v> = <u[0]e[0]+...+u[n-1]e[n-1], v[0]e[0]+...+u[n-1]e[n-1]>         //
//          = u[0]v[0] <e[0],e[0]> + ... + u[0]v[n-1] <e[n-1],e[n-1]>         //
//           + ... +                                                          //
//           u[n-1]v[0] <e[n-1,e[0]> + ... + u[n-1]v[n-1] <e[n-1],e[n-1]>     //
//          =  u[0]v[0] + ... + u[n-1]v[n-1]                                  //
//                                                                            //
//     The arguments u and v should be declared as double u[N] and            //
//     double v[N] where N >= n in the calling program.                       //
//                                                                            //
//  Arguments:                                                                //
//     double u[]  Pointer to the first element of the vector u.              //
//     double v[]  Pointer to the first element of the vector v.              //
//     int     n   The number of components of the vectors u and v.           //
//                                                                            //
//  Return Values:                                                            //
//     Inner Product of u and v.                                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double u[N], v[N], inner_product;                                      //
//                                                                            //
//     (your code to intialize the vectors u and v)                           //
//     inner_product = Inner_Product(u,v,N);                                  //
//     printf(" <u,v> = %12.6f\n", inner_product);                            //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Inner_Product(double u[], double v[], int n)
{
   double inner_product = 0.0;
   for (n--; n >= 0; n--) inner_product +=  u[n] * v[n];
   return inner_product;
}
////////////////////////////////////////////////////////////////////////////////
// File: vector_max_norm.c                                                    //
// Routines:                                                                  //
//    Vector_Max_Norm                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Vector_Max_Norm(double v[], int n)                                 //
//                                                                            //
//  Description:                                                              //
//     Given an n-d real vector space V with a basis consisting of vectors    //
//     with unit norm, the max norm on the vector space V is the maximum of   //
//     the absolute values of the components of a vector v with respect to    //
//     that basis i.e. for v = v[0]e[0] + ... + v[n-1]e[n-1], where e[0], ...,//
//     e[n-1] are the basis vectors for which                                 //
//     || e[0] || = ... = || e[n-1] || = 1, then                              //
//                   || v || = Max( |v[0]|, ..., |v[n-1]| ).                  //
//                                                                            //
//  Arguments:                                                                //
//     double v[]  Pointer to the first element of the vector v[n].           //
//     int     n   The number of elements of the vector v[].                  //
//                                                                            //
//  Return Values:                                                            //
//     max norm of v.                                                         //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N], norm;                                                     //
//                                                                            //
//     (your code to initialize the vector v)                                 //
//     norm = Vector_Max_Norm(v, N);                                          //
//     printf(" || v || = %12.6f\n", norm);                                   //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Vector_Max_Norm(double v[], int n)
{
   double norm = 0.0;
   double x;
   int i;
   for (i = 0; i < n; i++) if (norm < ( x = fabs( v[i] ) ) ) norm = x;
   return norm;
}
////////////////////////////////////////////////////////////////////////////////
// File: vector_l2_norm.c                                                     //
// Routines:                                                                  //
//    Vector_L2_Norm                                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Vector_L2_Norm(double v[], int n)                                  //
//                                                                            //
//  Description:                                                              //
//     Given an n-d real vector space V with a basis consisting of vectors    //
//     with unit norm, the l2 norm on the vector space V is the sqrt of the   //
//     sum of squares of the components of a vector v with respect to that    //
//     basis i.e.                                                             //
//     for v = v[0]e[0] + ... + v[n-1]e[n-1], where e[0], ..., e[n-1] are     //
//     the basis vectors for which || e[0] || = ... = || e[n-1] || = 1,       //
//     then                                                                   //
//               || v || = sqrt( |v[0]|^2  + ... + |v[n-1]|^2 ).              //
//                                                                            //
//  Arguments:                                                                //
//     double v[]  Pointer to the first element of the vector v[n].           //
//     int     n   The number of elements of the vector v[].                  //
//                                                                            //
//  Return Values:                                                            //
//     l2 norm of v.                                                          //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N], norm;                                                     //
//                                                                            //
//     (your code to intialize the vector v)                                  //
//     norm = Vector_L2_Norm(v, N);                                           //
//     printf(" || v || = %12.6f\n", norm);                                   //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Vector_L2_Norm(double v[], int n)
{
   double norm = 0.0;
   int i;
   for (i = 0; i < n; i++) norm +=  v[i] * v[i];
   return sqrt(norm);
}
////////////////////////////////////////////////////////////////////////////////
// File: div_vector_by_scalar.c                                               //
// Routine(s):                                                                //
//    Divide_Vector_by_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Divide_Vector_by_Scalar(double *v, double x, int n)                  //
//                                                                            //
//  Description:                                                              //
//     Divide the vector v by the non-zero scalar x, i.e. divide each         //
//     component of the vector v by the scalar x, v[i] <- v[i] / x for all i. //
//                                                                            //
//  Arguments:                                                                //
//     double *v    Pointer to the first element of the vector v.             //
//     double x     Scalar which divides each element of the vector v.        //
//     int    n     The number of components of the vector v.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N],  x;                                                       //
//                                                                            //
//     (your code to initialize the vector v and scalar x)                    //
//                                                                            //
//     if ( x != 0.0)  Divide_Vector_by_Scalar(v, x,N);                       //
//      printf("The vector v is \n"); ...                                     //
////////////////////////////////////////////////////////////////////////////////
void Divide_Vector_by_Scalar(double v[], double x, int n)
{
   double z = 1.0 / x;
   for (; n > 0; n--) *v++ *= z;
}

/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
void PrintVector2D(const vector<double> A[], int num_rows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(5);
  for (int i = 0; i < num_rows; i++) {
    for (
      auto it = A[i].begin();
        it != A[i].end(); it++)
        cout << *it << " ";
      cout << endl;
    }  
}

void PrintVector(const double TempVector[], int num_rows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(5);
  for (int i = 0; i < num_rows; i++)
    std::cout << TempVector[i] << "   "; 
  std::cout << endl;
}

int ReadMatrixMarketFile(matrix<double>& A, int& num_rows, int& num_cols) {

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

// ReadVectorMarketFile implementation for vectors
int ReadVectorMarketFile(vector<double>& B, int& num_cols) {
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

int Inverse_Power_Input()
{
	int i {0}, num_rows {0}, num_cols {0};
    matrix<double> A;
    ReadMatrixMarketFile(A,num_rows,num_cols);
    vector <double> BVector;
    ReadVectorMarketFile(BVector, num_cols);
    vector <double> array_A;  
    int max_tries {1000};
    double tolerance {0.00001};
    for (int i = 0; i < num_rows; i++)
    {
      for (int j = 0; j < num_cols; j++)
      {
	    array_A.push_back(A[i][j]);
	  }   
    }    
    vector <double> B;
    vector<double> X = CreateZeroVector(num_rows);
    for (i = 0; i < num_rows; i++)
    {
      B.push_back(BVector[i]);
	}
    // Place Parameter_Solution function call here.
    double eigenvalue {0.0};
    Inverse_Power_Solve (&array_A[0], num_cols, eigenvalue, &B[0], tolerance, max_tries); 
    array_A.clear();
    array_A.shrink_to_fit();
    BVector.clear();
    BVector.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    return 0;
}


int main ()
{
    Inverse_Power_Input();
    return 0;
}
