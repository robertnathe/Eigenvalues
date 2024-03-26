#include "Input_Template.h"
#include "Input_Template.cpp"

#include <cfloat>
using namespace std; 
int Rayleigh_Eigenvalue_Solve(vector<double> A[], double X[], double B[], unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter);
int Eigenvalue_Rayleigh_Method(vector<double> A[], unsigned int num_cols, double eigenvalue, double X[], double B[],double tolerance, unsigned int max_iter);
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols);
double Inner_Product(double u[], double v[], unsigned int num_rows);
double Vector_L2_Norm(double v[], unsigned int num_cols);
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int num_cols); 
int Output_Data(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter);
void PrintVector2D(const vector<double> A[], unsigned int num_rows);
void PrintVector(const double B[], unsigned int num_rows);
 
template<typename T>
using matrix = std::vector< std::vector<T> >;

int Rayleigh_Eigenvalue_Solve(vector<double> A[], double X[], double B[], unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter) {

  int err;
  err = Eigenvalue_Rayleigh_Method(&A[0], num_cols, eigenvalue, &X[0], &B[0], tolerance, max_iter);
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
int Eigenvalue_Rayleigh_Method(vector<double> A[], unsigned int num_cols, double eigenvalue, double X[], double B[],double tolerance, unsigned int max_iter)
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
   norm = Vector_L2_Norm(px0, num_cols);
   if ( norm == 0.0 )  return -2;
                       // Normalize the initial vector.
   Divide_Vector_by_Scalar(px0, norm, num_cols);
           // Calculate initial estimate of the dominant eigenvalue.
   Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], num_cols, num_cols);
   eigenvalue = Inner_Product(px, px0, num_cols);
      // Iterate for a maximum of max_iter times using the power method.
      // Exit if an error occurs or if the relative difference of two
      // successive estimates is less than the specified tolerance.
   for (tries = 1; tries < max_iter; tries++) {
      old_estimate = eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      norm = Vector_L2_Norm(px0, num_cols);
      if ( norm == 0.0 )  return -2;
      Divide_Vector_by_Scalar(px0, norm, num_cols);
      Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], num_cols, num_cols);
   eigenvalue = Inner_Product(px,px0,num_cols);
   if (eigenvalue == 0.0) return -3;
   if (fabs((eigenvalue - old_estimate) / eigenvalue) < tolerance)
   {
     revised_lambda = old_estimate;
     double norm_of_eigenvalue = Vector_L2_Norm(px,num_cols);
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
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols)
{
   unsigned int i,j;
   double sum = 0.0;
   for (i = 0; i < num_rows; i++)
   {
     sum = 0.0;
     for (j = 0; j < num_cols; j++)
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
double Inner_Product(double U[], double V[], unsigned int num_rows)
{
   double Inner_Product {0.0};
   for (num_rows--; num_rows > 0; num_rows--) {
	   Inner_Product +=  U[num_rows] * V[num_rows];
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
double Vector_L2_Norm(double v[], unsigned int num_cols)
{
   double norm = 0.0;
   unsigned int i;
   for (i = 0; i < num_cols; i++) norm +=  v[i] * v[i];
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
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int num_cols) 
{
   double z = 1.0 / TempScalar;
   for (; num_cols > 0; num_cols--) *X++ *= z;
}
int Output_Data(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter)
{
  unsigned int i;
  std::cout << std::setprecision (7) << endl;
  std::cout << "******************** Solve Ax = B ********************" << endl;
  // Displaying the 2D vector
  std::cout << "The vector A is the following: " << endl; 
  PrintVector2D(&A[0], num_cols);  
  std::cout << "The vector X is the following: " << endl;
  PrintVector(&X[0], num_cols);
  std::cout << "The vector B is the following: " << endl; 
  PrintVector(&B[0], num_cols); 
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << num_cols << " " << num_rows;
    outputFile << endl;
    for (i = 0; i < num_cols; i++)
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
void PrintVector2D(const vector<double> A[], unsigned int num_rows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(7);
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
  std::cout << std::setprecision(7);
  for (unsigned int i = 0; i < num_rows; i++)
    std::cout << TempVector[i] << "   "; 
  std::cout << endl;
}


int main ()
{
   //double omega {1.25};
   //(void) omega;
   int max_iter {1000};
   double eigenvalue {0.0};
   double tolerance {0.00001};
   std::cout << std::setprecision(5);     
   int num_rows {0}, num_cols {0};
   matrix<double> A;
   ReadMatrixMarketFile(A,num_rows,num_cols);
   vector <double> B;
   ReadVectorMarketFile(B, num_cols);      
   vector<double> X = CreateZeroVector(num_rows);
   vector<double> array_A = CreateZeroVector(num_rows*num_cols);
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
   
   ReadTridiagonal(A, B, X, num_rows, num_cols);
   for (int j = 0; j < num_cols; j++)
   {
     X[j] = 0.0; 
   }
   ReadVectorMarketFile(B, num_cols);  
   
   // Using time point and system_clock
   std::chrono::time_point<std::chrono::system_clock> start, end;
 
   start = std::chrono::system_clock::now();	 
   // Begin to place new code here.
   Rayleigh_Eigenvalue_Solve (&A[0],  &X[0], &B[0], num_cols, eigenvalue, tolerance, max_iter); 
   Output_Data(&A[0], &X[0], &B[0], num_rows, num_cols, eigenvalue, tolerance, max_iter);     
   
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
