#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()
using namespace std; 
#define a(i,j) a[(i)*nrows+(j)]
int Rayleigh_Eigenvalue_Solve(vector<double> A[], double X[], double B[], unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);
int Eigenvalue_Rayleigh_Method(vector<double> A[], unsigned int ncols, double eigenvalue, double X[], double B[],double tolerance, unsigned int max_iter);
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
double Inner_Product(double u[], double v[], unsigned int nrows);
double Vector_L2_Norm(double v[], unsigned int ncols);
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols); 
int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);
void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double B[], unsigned int nrows);
 
template<typename T>
using matrix = std::vector< std::vector<T> >;

int Rayleigh_Eigenvalue_Solve(vector<double> A[], double X[], double B[], unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter) {

  int err;
  err = Eigenvalue_Rayleigh_Method(&A[0], ncols, eigenvalue, &X[0], &B[0], tolerance, max_iter);
  if (err == -1) {
    printf(" The iteration failed to converge\n");
  }
  if (err == -2)
  {
    printf(" The iteration vector vanished\n");
  }
  if (err == -3)
  {
    printf(" The estimate of the dominant eigenvalue vanished\n");
  }

  return 0;	
}
int Eigenvalue_Rayleigh_Method(vector<double> A[], unsigned int ncols, double eigenvalue, double X[], double B[],double tolerance, unsigned int max_iter)
{
   unsigned int tries;
   double old_estimate = 0.0;
   double norm;
   double *px = X, *px0 = B, *pdum;
   double revised_lambda;
   (void) revised_lambda;
        // Check that the user specified tolerance is greater than the
        // machine epsilon.  If not, set it to the machine epsilon.
   if (tolerance < DBL_EPSILON) tolerance = DBL_EPSILON;
                // Check that the initial vector is non-zero.
   norm = Vector_L2_Norm(px0, ncols);
   if ( norm == 0.0 )  return -2;
                       // Normalize the initial vector.
   Divide_Vector_by_Scalar(px0, norm, ncols);
           // Calculate initial estimate of the dominant eigenvalue.
   Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], ncols, ncols);
   eigenvalue = Inner_Product(px, px0, ncols);
      // Iterate for a maximum of max_iter times using the power method.
      // Exit if an error occurs or if the relative difference of two
      // successive estimates is less than the specified tolerance.
   for (tries = 1; tries < max_iter; tries++) {
      old_estimate = eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      norm = Vector_L2_Norm(px0, ncols);
      if ( norm == 0.0 )  return -2;
      Divide_Vector_by_Scalar(px0, norm, ncols);
      Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], ncols, ncols);
   eigenvalue = Inner_Product(px,px0,ncols);
   if (eigenvalue == 0.0) return -3;
   if (fabs((eigenvalue - old_estimate) / eigenvalue) < tolerance)
   {
     revised_lambda = old_estimate;
     double norm_of_eigenvalue = Vector_L2_Norm(px,ncols);
     std::cout << "The norm of the eigenvalue is the following: " << norm_of_eigenvalue << endl;
     break;
   }
  }
  if (tries >= max_iter) return -1;
  return tries;
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
   double Inner_Product {0.0};
   for (nrows--; nrows > 0; nrows--) {
	   Inner_Product +=  U[nrows] * V[nrows];
   }
   return Inner_Product;
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
double Vector_L2_Norm(double v[], unsigned int ncols)
{
   double norm = 0.0;
   unsigned int i;
   for (i = 0; i < ncols; i++) norm +=  v[i] * v[i];
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
    double value, elem {1.0};
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
      X.push_back(1.0);
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
    Rayleigh_Eigenvalue_Solve (&A[0],  &X[0], &B[0], ncols, eigenvalue, tolerance, max_iter);
//    print_results (&array_A[0], ncols, &B[0], tolerance, max_iter);
//    Power_Eigenvalue_Solve (&A[0], &X[0], &B[0], ncols, eigenvalue, tolerance, max_iter);  
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
