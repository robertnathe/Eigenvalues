#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()
using namespace std; 
#define a(i,j) a[(i)*nrows+(j)]
int Power_Eigenvalue_Solve (vector<double> A[], double X[], double B[], unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);
int Eigenvalue_Power_Method (vector<double> A[], unsigned int nrows, double eigenvalue, double X[], double X0[], double tolerance, unsigned int max_iter);
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
double Inner_Product(double u[], double v[], unsigned int nrows);
double Vector_Max_Norm(double TempVector[], unsigned int ncols);
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols); 
int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);
void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double B[], unsigned int nrows);
  
template<typename T>
using matrix = std::vector< std::vector<T> >;

int Power_Eigenvalue_Solve (vector<double> A[], double X[], double B[], unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter)
{
  //by default all values are 0
  // Enter code to process data here.
  unsigned int tries = { 1 };
  int err = {0};
  std::cout << "Prog: Power_Method_for_Eigenvalues.cpp" << endl << endl << endl;
  err = Eigenvalue_Power_Method (&A[0], ncols, eigenvalue, &X[0], &B[0], tolerance, max_iter);
  if (tries >= max_iter)
    return -1;
  if (err == -1)
    std::cout << " The iteration failed to converge" << endl;
  else if (err == -2)
    std::cout << " The iteration vector vanished" << endl;
  else if (err == -3)
    std::cout << " The estimate of the dominant eigenvalue vanished" << endl;
  return 0;
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
int Eigenvalue_Power_Method (vector<double> A[], unsigned int nrows, double eigenvalue, double X[],
			 double X0[], double tolerance, unsigned int max_iter)
{
  unsigned int tries;
  double old_estimate;
  double norm;
  double *px = X, *px0 = X0, *pdum;
  // Check that the user specified tolerance is greater than the
  // machine epsilon.  If not, set it to the machine epsilon.
  if (tolerance < DBL_EPSILON)
    tolerance = DBL_EPSILON;
  // Check that the initial vector is non-zero.
  // norm = Vector_L2_Norm(px0, nrows);
  norm = Vector_Max_Norm (px0, nrows);
  if (norm == 0.0)
    return -2;
  // Normalize the initial vector.
  Divide_Vector_by_Scalar (px0, norm, nrows);
  // Calculate initial estimate of the dominant eigenvalue.
  Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], nrows, nrows);
  eigenvalue = Inner_Product (px, px0, nrows);
  // Iterate for a maximum of max_tries times using the power method.
  // Exit if an error occurs or if the relative difference of two
  // successive estimates is less than the specified tolerance.
  for (tries = 1; tries < max_iter; tries++)
    {
      old_estimate = eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      norm = Vector_Max_Norm (px0, nrows);
      if (norm == 0.0)
	return -2;
      Divide_Vector_by_Scalar (px0, norm, nrows);
      Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], nrows, nrows);
      eigenvalue = Inner_Product (px, px0, nrows);
      if (eigenvalue == 0.0)
	return -3;
      if (fabs ((eigenvalue - old_estimate) / eigenvalue) < tolerance)
       {
		 std::cout << "The norm of the dominant eigenvalue is the following: " << norm << endl;
         break;
       }
    }
  if (tries >= max_iter)
  {
    return -1;
  }
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
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   double sum = 0.0;
   for (i = 0; i < nrows; i++)
   {
     sum = 0.0;
     for (j = 0; j < ncols; j++)
     {
       sum += A[i][j] * X[j];
       B[i] = sum;
     }
   }
   return 0;
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
double Inner_Product(double U[], double V[], unsigned int nrows)
{
   double Inner_Product = 0.0;
   for (nrows--; nrows > 0; nrows--) Inner_Product +=  U[nrows] * V[nrows];
   return Inner_Product;
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
double Vector_Max_Norm(double TempVector[], unsigned int ncols)
{
   double norm = 0.0;
   double TempScalar;
   unsigned int i;
   for (i = 0; i < ncols; i++) if (norm < ( TempScalar = fabs( TempVector[i] ) ) ) norm = TempScalar;
   return norm;
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
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols) 
{
   double z = 1.0 / TempScalar;
   for (; ncols > 0; ncols--) *X++ *= z;
}
int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter)
{
  unsigned int i;
  std::cout << std::setprecision (7) << endl;
  std::cout << "******************** Solve Ax = B ********************" << endl;
  // Displaying the 2D vector
  std::cout << "The vector A is the following: " << endl; 
  PrintVector2D(&A[0], ncols);  
  std::cout << "The vector X is the following: " << endl;
  PrintVector(&X[0], ncols);
  std::cout << "The vector B is the following: " << endl; 
  PrintVector(&B[0], ncols); 
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << ncols << " " << nrows;
    outputFile << endl;
    for (i = 0; i < ncols; i++)
      outputFile <<  1 << " " << i+1 << " " << X[i] << endl;
    // Close the file
    outputFile.close();
  } else {
    std::cout << "Error writing to file." << std::endl;
  }
  return 0;
}
/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
void PrintVector2D(const vector<double> A[], unsigned int nrows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(7);
  for (unsigned int i = 0; i < nrows; i++) {
    for (
      auto it = A[i].begin();
        it != A[i].end(); it++)
        cout << *it << " ";
      cout << endl;
    }  
}
void PrintVector(const double TempVector[], unsigned int nrows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(7);
  for (unsigned int i = 0; i < nrows; i++)
    std::cout << TempVector[i] << "   "; 
  std::cout << endl;
}
class Input {
  public: 
  unsigned int Input_Begin(void)
  {
    unsigned int i,j, max_iter {1000};
    double tolerance {0.00001};
    double eigenvalue {0.0};
    FILE *myFile;
    myFile = fopen ("A.dat", "r");
    unsigned int dimArray[3];
    unsigned int nrows, number_of_entries_A;
    unsigned int ncols;
    unsigned int i_index, j_index;
    double value, elem {0.0};
    while (myFile == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }
    fscanf (myFile, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
    fscanf (myFile, "%u,", &dimArray[i]);
    }
    nrows = dimArray[0];
    ncols = dimArray[1];
    number_of_entries_A = dimArray[2];
    vector <double> array_A;
    for (i = 0; i < number_of_entries_A; i++)
    {
      fscanf (myFile, "%u,", &i_index);
      i_index--;
      fscanf (myFile, "%u,", &j_index);
      j_index--;
      fscanf (myFile, "%lf,", &value);
//    Change program to use the single index array_A
      array_A.push_back(value);
    }
    fclose (myFile);
    FILE *myFile2;
    myFile2 = fopen ("B.dat", "r");
    unsigned int dim_B_Array[3];
    unsigned int number_of_entries_B;
    unsigned int col_B;
    (void) col_B;
    while (myFile2 == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }  
    fscanf (myFile2, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
      fscanf (myFile2, "%u,", &dim_B_Array[i]);
    }
    col_B = dim_B_Array[1];
    number_of_entries_B = dim_B_Array[2];
    vector <double> B;
    vector <double> X;
    for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%u,", &i_index);
      i_index--;
      fscanf (myFile2, "%u,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B.push_back(value);
      X.push_back(0.0);
    }
    fclose (myFile2);
    // Initializing the vector of vectors
//    vector<vector<double> > A;
    matrix<double> A;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> v1;  
        for (j = 0; j < ncols; j++) {
		  elem = array_A[i*nrows+j];
          v1.push_back(elem);           
            if (j == (ncols)) {
				v1.push_back(elem);
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A.push_back(v1);
    }    
    Power_Eigenvalue_Solve (&A[0], &X[0], &B[0], ncols, eigenvalue, tolerance, max_iter);  
    Output_Data(&A[0], &X[0], &B[0], nrows, ncols, eigenvalue, tolerance, max_iter);    
    array_A.clear();
    array_A.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    A.clear();
    A.shrink_to_fit();
    return 0;
  }  
};

int main ()
{
  Input Input_One;
  Input_One.Input_Begin();
  return 0;
}
