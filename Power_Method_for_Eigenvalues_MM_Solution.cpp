#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()
using namespace std;
int Eigenvalue_Power_Method (double *A, int n, double *eigenvalue, double x[],
			     double x0[], double tolerance, int max_tries);
void Divide_Vector_by_Scalar (double v[], double x, int n);
double Vector_Max_Norm (double v[], int n);
double Inner_Product (double u[], double v[], int n);
void Multiply_Matrix_by_Vector (double u[], double *A, int nrows, int ncols,
				double v[]);
int Input_matrix_A (void);
int Parameters_Solve (double array_A[], double B[], double X[], int row,
		      int col);
int print_results (double array[], double x0[], double x[], int row, int col);
int
main ()
{
  Input_matrix_A ();
  return 0;
}

int
Input_matrix_A (void)
{
  FILE *myFile;
  myFile = fopen ("A.dat", "r");
  //read file into array
  int dimArray[3];
  int row, col, number_of_entries;
  int i_index, j_index;
  double value;
  int i;
  if (myFile == NULL)
    {
      printf ("Error Reading File\n");
      exit (0);
    }
  fscanf (myFile, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
    {
      fscanf (myFile, "%d,", &dimArray[i]);
    }
  row = dimArray[0];
  col = dimArray[1];
  number_of_entries = dimArray[2];
  double A[row][col];
//    print_matrix( (double*) A, row, col);
  for (i = 0; i < number_of_entries; i++)
    {

      fscanf (myFile, "%d,", &i_index);
      i_index--;
      fscanf (myFile, "%d,", &j_index);
      j_index--;
      fscanf (myFile, "%lf,", &value);
      A[i_index][j_index] = value;
    }
  fclose (myFile);
  FILE *myFile2;
  myFile2 = fopen ("B.dat", "r");
  //read file into array
  int dim_B_Array[3];
  int row_B, col_B, number_of_entries_B;
  if (myFile2 == NULL)
    {
      printf ("Error Reading File\n");
      exit (0);
    }
  fscanf (myFile2, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
    {
      fscanf (myFile2, "%d,", &dim_B_Array[i]);
    }
  row_B = dim_B_Array[0];
  col_B = dim_B_Array[1];
  number_of_entries_B = dim_B_Array[2];
  double B_temporary[row_B][col_B];
  double B[col_B], X[col_B];
  for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%d,", &i_index);
      i_index--;
      fscanf (myFile2, "%d,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B_temporary[i_index][j_index] = value;
      B[j_index] = B_temporary[i_index][j_index];
    }
  fclose (myFile2);
  int j, kk;
  kk = row * col;
  double array_A[kk];
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  array_A[i * row + j] = A[i][j];
	}
    }
  printf ("Prog: Power_MM_Solution.cpp\n");
  Parameters_Solve (&array_A[0], &B[0], &X[0], row, col);
  return 0;
}

int
Parameters_Solve (double array_A[], double B[], double X[], int row, int col)
{
  int i, j;
  double A[row][col];
//  int pivot[col] = { 0 };
  B[col] =
  {
  0};
  X[col] =
  {
  0};
  fprintf (stdout, "\n");
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  A[i][j] = array_A[i * row + j];
	}
    }
  // Enter code to process data here.
  double x[col] = { 0.0, 0.0, 0.0 };
  double x0[col] = { 1.0, 1.0, 1.0 };
  double eigenvalue = { 0.0 };
  double tolerance = 0.000001;
  int max_tries = 30;
  int tries = { 1 }, err =
  {
  0};
  printf ("Prog: Power_Solution.cpp\n\n\n");
//  print_results((double*)A, &x[0], &x0[0]);
  Eigenvalue_Power_Method ((double *) A, col, &eigenvalue, &x[0], &x0[0],
			   tolerance, max_tries);
  if (tries >= max_tries)
    return -1;
  err =
    Eigenvalue_Power_Method ((double *) A, col, &eigenvalue, x, x0, tolerance,
			     max_tries);
  if (err == -1)
    printf (" The iteration failed to converge\n");
  else if (err == -2)
    printf (" The iteration vector vanished\n");
  else if (err == -3)
    printf (" The estimate of the dominant eigenvalue vanished\n");
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  array_A[i * row + j] = A[i][j];
	}
    }
  print_results (&array_A[0], &x0[0], &x[0], row, col);
  return 0;
}

int
print_results (double array[], double x0[], double x[], int row, int col)
{
  int i, j;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  printf ("%6.3f   ", array[i * row + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("and x0 = \n");
  for (i = 0; i < col; i++)
    printf ("%6.3f   ", x0[i]);
  printf ("\n\n");
  printf ("The solution is x = \n");
  for (i = 0; i < col; i++)
    {
      printf ("%6.3f   ", x[i]);
      // Assigning values to the arrays X and pivot to fix warning message
    }
  printf ("\n\n");
  printf ("\n\n\n");
  return 0;
}				// static void

// print_matrix (double *X, int m, int n)
// {
//  int i, j;
//  for (i = 0; i < m; i++)
//    {
//      for (j = 0; j < n; j++)
//      printf (" %8.2f ", *X++);
//      printf ("\n");
//    }
// }

// static void print_vector (double *B, int n)
// {
//  printf ("\nThe  matrix column B is the following: \n\n");
//  int i;
//  for (i = 0; i < n; i++)
//    printf (" %8.2f ", B[i]);
//  printf ("\n");
// }

////////////////////////////////////////////////////////////////////////////////
// File: zero_matrix.c                                                        //
// Routine(s):                                                                //
//    Zero_Matrix                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Zero_Matrix(double *A, int nrows, int ncols)                         //
//                                                                            //
//  Description:                                                              //
//     Set the nrows x ncols matrix A equal to the zero matrix, i.e.          //
//     A[i][j] = 0 for all i, j.                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A.                       //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[N][M];                                                        //
//                                                                            //
//     Zero_Matrix(&A[0][0], N, M);                                           //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
// int Zero_Matrix (double *X, int nrows, int ncols)
// {
//  int n = nrows * ncols;
//  for (; n > 0; n--)
//    *X++ = 0.0;
//  return 0;
// }
////////////////////////////////////////////////////////////////////////////////
// File: zero_vector.c                                                        //
// Routine(s):                                                                //
//    Zero_Vector                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Zero_Vector(double *A, int n)                                        //
//                                                                            //
//  Description:                                                              //
//     Set the vector A equal to the zero vector, i.e. A[i] = 0 for all i.    //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the vector A.             //
//     int    n     The number of components of the vector A.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N];                                                           //
//                                                                            //
//     Zero_Vector(A, N);                                                     //
//     printf("The vector A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
// void Zero_Vector (double *X, int n)
// {
//  for (; n > 0; n--)
//    *X++ = 0.0;
// }
int
Eigenvalue_Power_Method (double *A, int n, double *eigenvalue, double x[],
			 double x0[], double tolerance, int max_tries)
{
  int tries;
  double old_estimate;
  double norm;
  double *px = x, *px0 = x0, *pdum;
  // Check that the user specified tolerance is greater than the 
  // machine epsilon.  If not, set it to the machine epsilon.
  if (tolerance < DBL_EPSILON)
    tolerance = DBL_EPSILON;
  // Check that the initial vector is non-zero.
  // norm = Vector_L2_Norm(px0, n);
  norm = Vector_Max_Norm (px0, n);
  if (norm == 0.0)
    return -2;
  // Normalize the initial vector.
  Divide_Vector_by_Scalar (px0, norm, n);
  // Calculate initial estimate of the dominant eigenvalue.
  Multiply_Matrix_by_Vector (px, A, n, n, px0);
  *eigenvalue = Inner_Product (px, px0, n);
  // Iterate for a maximum of max_tries times using the power method.
  // Exit if an error occurs or if the relative difference of two
  // successive estimates is less than the specified tolerance.
  for (tries = 1; tries < max_tries; tries++)
    {
      old_estimate = *eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      norm = Vector_Max_Norm (px0, n);
      double dominant_eigenvalue = *eigenvalue;
      if (norm == 0.0)
	return -2;
      Divide_Vector_by_Scalar (px0, norm, n);
      Multiply_Matrix_by_Vector (px, A, n, n, px0);
      *eigenvalue = Inner_Product (px, px0, n);
      if (*eigenvalue == 0.0)
	return -3;
      if (fabs ((*eigenvalue - old_estimate) / *eigenvalue) < tolerance)
	break;
    }
  if (tries >= max_tries)
    return -1;
    printf("The norm of the dominant eigenvalue is %lf.\n", norm);
  return tries;
}

////////////////////////////////////////////////////////////////////////////////
//  int Eigenvalue_Power_Method(double *A, int n, double *eigenvalue,         //
//                 double x[], double x0[], double tolerance, int max_tries)  //
//                                                                            //
//  Description:                                                              //
//     Let A be a real n x n matrix with n linearly independent eigenvectors, //
//     the power method is a method for estimating the dominant eigenvalue of //
//     A.  In particular, given an arbitrary starting vector x0, x0 can be    //
//     written as x0 = a[0] u[0] + ... + a[n-1] u[n-1], where {u[i]} are the  //
//     n linearly independent eigenvectors and a[i] are scalars.  Applying    //
//     A to x0, A x0 = a[0] z[0] u[0] + ... + a[n-1] z[n-1] u[n-1],  where    //
//     z[i] is the eigenvalue with corresponding eigenvector u[i].  Let x[k]  //
//     be the result of applying A to x0 k times.  Then x[k] = A x[k-1], and  //
//     if z[0] is the dominant eigenvalue, as k becomes large, x[k] approaches//
//     z[0]^k {a[0] u[0] + ... + a[n-1] (z[n-1]/z[0])^k u[n-1]} which is      //
//     approximately z[0]^k a[0] u[0] ~ z[0] x[k-1].  Since |z[0]|^k either   //
//     diverges as k becomes large or tends to 0 as k becomes large, at each  //
//     step if x[k] is normalized, then z[0] ~ x[k]'x[k-1] / x[k-1]'x[k-1].   //
//                                                                            //
//     The power method is:  (1) Start with an initial non-zero vector x0,    //
//     (2) Normalize the current vector x[k], (3) x[k+1] = A x[k],            //
//     (4) let z = x[k+1]'x[k], (5) if the relative difference of the new     //
//     value of z and the old value of z is less than a preassigned tolerance,//
//     then halt the procedure; otherwise go to step (2).                     //
//                                                                            //
//     Note! If the matrix A does not admit n linearly independent eigen-     //
//     vectors, a dominant eigenvalue may still be found.                     //
//                                                                            //
//     Note! If the dominant eigenvalue is not unique, the process may not    //
//     converge.  If z and -z are dominant eigenvalues, then the process      //
//     defined by x[k+2]'x[k] / x[k]'x[k] converges to z^2.  If the dominant  //
//     eigenvalue is complex, then the complex conjugate is also a dominant   //
//     eigenvalue.  And some refinements to this procedure need to be made.   //
//                                                                            //
//     Note! If the initial vector x0 belongs to the orthogonal complement    //
//     of the eigenvector corresponding to the dominant eigenvalue, then      //
//     theoretically all subsequent iterates will also.  But round-off errors //
//     will contaminate this theoretical result and introduce non-zero        //
//     components of x[k] along the direction of the eigenvector corresponding//
//     to the dominant eigenvalue.  Convergence then behaves as before.       //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *A                                                              //
//        The pointer to the first element of the matrix A[n][n].  The matrix //
//        A is not altered during the computations.                           //
//     int     n                                                              //
//        The number of rows and/or columns of the matrix A.                  //
//     double *eigenvalue                                                     //
//        If the call returns a positive number indicating success, then the  //
//        estimate of the eigenvalue is set in the address pointed to by      //
//        *eigenvalue.                                                        //
//     double x[]                                                             //
//        Working storage, should be declared in the calling routine as       //
//        "double x[n]", where n is as above.                                 //
//     double x0[]                                                            //
//        On input the initial vector.  It should be declared in the calling  //
//        routine as "double x0[n]", where n is as above.  It should be       //
//        initialized as well.                                                //
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
//      -3    Failure - Estimate of dominant eigenvalue is 0.                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], x[N], x0[N];                                           //
//     double eigenvalue, tolerance;                                          //
//     int max_tries;                                                         //
//                                                                            //
//     (your code to initialize the matrix A, the vector x0, tolerance,       //
//      and max_tries)                                                        //
//     err = Eigenvalue_Power_Method((double*)A, N, &eigenvalue, x, x0,       //
//                                                     tolerance, max_tries); //
//     if (err == -1)                                                         //
//        printf(" The iteration failed to converge\n");                      //
//     else if (err == -2)                                                    //
//        printf(" The iteration vector vanished\n");                         //
//     else if (err == -3)                                                    //
//        printf(" The estimate of the dominant eigenvalue vanished\n");      //
//     else printf(" The eigenvalue is \n");                                  //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
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
void
Divide_Vector_by_Scalar (double v[], double x, int n)
{
  double z = 1.0 / x;
  for (; n > 0; n--)
    *v++ *= z;
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
double
Vector_Max_Norm (double v[], int n)
{
  double norm = 0.0;
  double x;
  int i;
  for (i = 0; i < n; i++)
    if (norm < (x = fabs (v[i])))
      norm = x;
  return norm;
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
double
Inner_Product (double u[], double v[], int n)
{
  double inner_product = 0.0;
  for (n--; n >= 0; n--)
    inner_product += u[n] * v[n];
  return inner_product;
}

////////////////////////////////////////////////////////////////////////////////
// File: multiply_matrix_by_vector.c                                          //
// Routine(s):                                                                //
//    Multiply_Matrix_by_Vector                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Matrix_by_Vector(double u[], double *A, int nrows,          //
//                                                   int ncols, double v[])   //
//                                                                            //
//  Description:                                                              //
//     Post multiply the nrows x ncols matrix A by the column vector v        //
//     to form the column vector u, i.e. u = A v.                             //
//     The matrix A should be declared as "double A[nrows][ncols]" in the     //
//     calling routine.  The vector v declared as "double v[ncols]" and       //
//     the vector u declared as "double u[nrows]" in the calling routine.     //
//                                                                            //
//  Arguments:                                                                //
//     double *u    Pointer to the first element of the vector u.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A and the number of      //
//                  components of the column vector u.                        //
//     int    ncols The number of columns of the matrices A and the           //
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
void
Multiply_Matrix_by_Vector (double u[], double *A, int nrows, int ncols,
			   double v[])
{
  int i, j;
  for (i = 0; i < nrows; A += ncols, i++)
    for (u[i] = 0.0, j = 0; j < ncols; j++)
      u[i] += A[j] * v[j];
}

////////////////////////////////////////////////////////////////////////////////
// File: eigen_power_method.c                                                 //
// Contents:                                                                  //
//    Eigenvalue_Power_Method                                                 //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    Multiply_Matrix_by_Vector                                               //
//    Inner_Product                                                           //
//    Vector_L2_Norm                                                          //
//    Divide_Vector_by_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
#include <float.h>		// required for DBL_EPSILON
#include <math.h>		// required for fabs() and for sqrt
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
