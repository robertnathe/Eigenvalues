#include<iomanip>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>		// required for fabs() and for sqrt()
#include <cmath>
#include <float.h>                  // required for DBL_EPSILON
using namespace std;
int Input_Data (void);
int Parameter_Solution(double array_A[], int ncols, double eigenvalue, double X[], double X0[],double tolerance, int max_tries);
int print_results (double array_A[], int ncols, double X0[],double tolerance, int max_tries);
//                    Required Externally Defined Routines
// Place Externally Define Routines here.
int array_Multiply_Matrix_by_Vector(double array_A[], double X[], double B[], int nrows, int ncols);
int Eigenvalue_Rayleigh_Method(double array_A[], int n, double eigenvalue, double X[], double X0[],double tolerance, int max_tries);
double Inner_Product(double u[], double v[], int n);
double Vector_L2_Norm(double v[], int n);
void   Divide_Vector_by_Scalar(double v[], double x, int n);

int main ()
{
  Input_Data();
  return 0;
}

int Eigenvalue_Rayleigh_Method(double array_A[], int n, double eigenvalue, double X[], double X0[],double tolerance, int max_tries)
{
   int tries;
   int nrows = n;
   int ncols = n;
   double old_estimate = 0.0;
   double norm;
   double *px = X, *px0 = X0, *pdum;
   double revised_lambda;
        // Check that the user specified tolerance is greater than the
        // machine epsilon.  If not, set it to the machine epsilon.
   if (tolerance < DBL_EPSILON) tolerance = DBL_EPSILON;
                // Check that the initial vector is non-zero.
   norm = Vector_L2_Norm(px0, n);
   if ( norm == 0.0 )  return -2;
                       // Normalize the initial vector.
   Divide_Vector_by_Scalar(px0, norm, n);
           // Calculate initial estimate of the dominant eigenvalue.
   array_Multiply_Matrix_by_Vector(&array_A[0], &X0[0], &X[0], nrows, ncols);
   eigenvalue = Inner_Product(px, px0, n);
      // Iterate for a maximum of max_tries times using the power method.
      // Exit if an error occurs or if the relative difference of two
      // successive estimates is less than the specified tolerance.
   for (tries = 1; tries < max_tries; tries++) {
      old_estimate = eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      norm = Vector_L2_Norm(px0, n);
      if ( norm == 0.0 )  return -2;
      Divide_Vector_by_Scalar(px0, norm, n);
      array_Multiply_Matrix_by_Vector(&array_A[0], &px0[0], &px[0], nrows, ncols);
   eigenvalue = Inner_Product(px,px0,n);
   if (eigenvalue == 0.0) return -3;
   if (fabs((eigenvalue - old_estimate) / eigenvalue) < tolerance)
   {
     revised_lambda = old_estimate;
     double dominant_eigenvalue = revised_lambda;
     printf("The dominant eigenvalue is %6.4lf .\n", dominant_eigenvalue);
     break;
   }
  }
  if (tries >= max_tries) return -1;
  return tries;
}

int array_Multiply_Matrix_by_Vector(double array_A[], double X[], double B[], int nrows, int ncols)
{
   int i,j;
   double sum = 0.0;
   for (i = 0; i < nrows; i++)
   {
     sum = 0.0;
     for (j = 0; j < ncols; j++)
     {
       sum += array_A[i*ncols+j] * X[j];
       B[i] = sum;
     }
   }
   return 0;
}

int Parameter_Solution(double array_A[], int ncols, double eigenvalue, double X[], double X0[],double tolerance, int max_tries)
{
  int err;
  err = Eigenvalue_Rayleigh_Method(&array_A[0], ncols, eigenvalue, &X[0], &X0[0], tolerance, max_tries);
  if (err == -1)
    printf(" The iteration failed to converge\n");
  if (err == -2)
  {
    printf(" The iteration vector vanished\n");
  }
  if (err == -3)
  {
    printf(" The estimate of the dominant eigenvalue vanished\n");
  }
  print_results (&array_A[0], ncols, &X0[0], tolerance, max_tries);
  return 0;
}

int Input_Data(void)
{
  int i, max_tries = 100;
  double eigenvalue {0.0};
  double tolerance {0.0001};
  FILE *myFile;
  myFile = fopen ("A.dat", "r");
  int dimArray[3];
  int nrows, ncols, number_of_entries_A;
  int i_index, j_index;
  double value;
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
  nrows = dimArray[0];
  ncols = dimArray[1];
  number_of_entries_A = dimArray[2];
  double *array_A;
  array_A = new double[nrows*ncols];
  double *X0;
  X0 = new double[ncols];
  for (i = 0; i < number_of_entries_A; i++)
  {
    fscanf (myFile, "%d,", &i_index);
    i_index--;
    fscanf (myFile, "%d,", &j_index);
    j_index--;
    fscanf (myFile, "%lf,", &value);
    array_A[i_index * ncols + j_index] = value;
  }
  fclose (myFile);
  FILE *myFile2;
  myFile2 = fopen ("B.dat", "r");
  int dim_X0_Array[3];
  int col_B, number_of_entries_X0;
  if (myFile2 == NULL)
  {
    printf ("Error Reading File\n");
    exit (0);
  }
  fscanf (myFile2, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
  {
    fscanf (myFile2, "%d,", &dim_X0_Array[i]);
  }
  col_B = dim_X0_Array[1];
  number_of_entries_X0 = dim_X0_Array[2];
  double *B;
  B = new double[col_B];
  double *X;
  X = new double[col_B];
  for (i = 0; i < number_of_entries_X0; i++)
  {
    fscanf (myFile2, "%d,", &i_index);
    i_index--;
    fscanf (myFile2, "%d,", &j_index);
    j_index--;
    fscanf (myFile2, "%lf,", &value);
    X0[j_index] = value;
  }
  fclose (myFile2);
  Parameter_Solution (&array_A[0], ncols, eigenvalue, &X[0], &X0[0], tolerance, max_tries);
  delete [] array_A;
  delete [] B;
  delete [] X;
  delete [] X0;
  return 0;
}

int print_results (double array_A[], int ncols, double X0[],double tolerance, int max_tries)
{
  int i, j;
  int nrows = ncols;
  printf ("The tolerance is ");
  printf ("%6.4lf \n", tolerance);
  printf ("The max_tries is %d \n", max_tries);
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where the matrix A = \n");
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
	{
	  printf ("%6.4f   ", array_A[i * ncols + j]);
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
  fprintf (myFile3,"%d %d %d\n", 1, nrows, ncols);
  printf ("The normalized eigenvector is X0 = \n");
  for (i = 0; i < ncols; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, X0[i]);
    printf ("%6.4f    ", X0[i]);
  }
  fclose(myFile3);
  printf ("\n\n");
  return 0;
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
void Multiply_Matrix_by_Vector(double u[], double *A, int nrows, int ncols, double v[])
{
   int i,j;
   for (i = 0; i < nrows; A += ncols, i++)
      for (u[i] = 0.0, j = 0; j < ncols; j++) u[i] += A[j] * v[j];
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
////////////////////////////////////////////////////////////////////////////////
// File: eigen_rayleigh_method.c                                              //
// Contents:                                                                  //
//    Eigenvalue_Rayleigh_Method                                              //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    Multiply_Matrix_by_Vector                                               //
//    Inner_Product                                                           //
//    Vector_L2_Norm                                                          //
//    Divide_Vector_by_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Eigenvalue_Rayleigh_Method(double *A, int n, double *eigenvalue,      //
//                 double x[], double x0[], double tolerance, int max_tries)  //
//                                                                            //
//  Description:                                                              //
//     Rayleigh's method is a variant of the power method in which the matrix //
//     A is a symmetric n x n real matrix.  A symmetric square n x n matrix   //
//     has n linearly independent eigenvectors, which may be chosen to be     //
//     orthonormal.  Similar to the power method Rayleigh's method is a method//
//     for estimating the dominant eigenvalue of the matrix A.                //
//     In particular, given an arbitrary starting vector x0, x0 can be written//
//     as  x0 = a[0] u[0] + ... + a[n-1] u[n-1], where {u[i]} are the n       //
//     orthonormal eigenvectors and a[i] are scalars.  Applying A to x0,      //
//            A x0 = a[0] z[0] u[0] + ... + a[n-1] z[n-1] u[n-1],             //
//     where z[i] is the eigenvalue with corresponding eigenvector u[i].  Let //
//     x[k] be the result of applying A to x0 k times.  Then x[k] = A x[k-1], //
//     and if z[0] is the dominant eigenvalue, as k becomes large, x[k]       //
//     approaches z[0]^k {a[0] u[0] + ... + a[n-1] (z[n-1]/z[0])^k u[n-1]}    //
//     which is approximately z[0]^k a[0] u[0] ~ z[0] x[k-1].  Since |z[0]|^k //
//     either diverges as k becomes large or tends to 0 as k becomes large, at//
//     each step if x[k-1] is normalized, then z[0] ~ x[k]'x[k-1].            //
//                                                                            //
//     Rayleigh's method is:  (1) Start with an initial non-zero vector x0,   //
//     (2) Normalize the current vector x[k], (3) x[k+1] = A x[k],            //
//     (4) let z = x[k+1]'x[k], (5) if the relative difference of the new     //
//     value of z and the old value of z is less than a preassigned tolerance,//
//     then halt the procedure; otherwise go to step (2).                     //
//                                                                            //
//     Note! If the dominant eigenvalue is not unique, then the process may   //
//     not converge.  If z and -z are dominant eigenvalues, then the process  //
//     defined by x[k+2]'x[k] / x[k]'x[k] converges to z^2.                   //
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
//     err = Eigenvalue_Rayleigh_Method((double*)A, N, &eigenvalue, x, x0,    //
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
