#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <chrono>
template<typename T>
using matrix = std::vector< std::vector<T> >;
using namespace std;
 
int Rayleigh_Eigenvalue_Solve(vector<double> A[], double X[], double B[], unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter);
int Eigenvalue_Rayleigh_Method(vector<double> A[], unsigned int num_cols, double eigenvalue, double X[], double B[],double tolerance, unsigned int max_iter);
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols);
double CalculateInnerProduct(double u[], double v[], unsigned int num_rows);
double Vector_L2_Norm(double v[], unsigned int num_cols);
void ScaleVector(double X[], double TempScalar, unsigned int num_cols); 
std::vector<double> CreateVectorFilledWithValue(int num_rows); 

int ReadMatrixMarketMatrix(matrix<double>& A, int& num_rows, int& num_cols);
int ReadMatrixMarketVector(vector<double>& B, int& num_cols); 
int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols);
int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols); 
int Write1DArrayToMatrixMarketFile(const double B[], int num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols);
void PrintMatrix(const vector<double> A[], unsigned int num_rows);
void PrintVector(const double B[], unsigned int num_rows);
void PrintVector(const std::vector<double>& vector);
void PrintMatrix(const matrix<double>& A); 
int Output_Data(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter);

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
   ScaleVector(px0, norm, num_cols);
           // Calculate initial estimate of the dominant eigenvalue.
   Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], num_cols, num_cols);
   eigenvalue = CalculateInnerProduct(px, px0, num_cols);
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
      ScaleVector(px0, norm, num_cols);
      Multiply_Matrix_by_Vector(&A[0], &px0[0], &px[0], num_cols, num_cols);
   eigenvalue = CalculateInnerProduct(px,px0,num_cols);
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
// File: CalculateInnerProduct.c                                                      //
// Routines:                                                                  //
//    CalculateInnerProduct                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double CalculateInnerProduct(double u[], double v[], int n)                       //
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
//     double u[N], v[N], CalculateInnerProduct;                                      //
//                                                                            //
//     (your code to intialize the vectors u and v)                           //
//     CalculateInnerProduct = CalculateInnerProduct(u,v,N);                                  //
//     printf(" <u,v> = %12.6f\n", CalculateInnerProduct);                            //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double CalculateInnerProduct(double U[], double V[], unsigned int num_rows)
{
   double CalculateInnerProduct {0.0};
   for (num_rows--; num_rows > 0; num_rows--) {
	   CalculateInnerProduct +=  U[num_rows] * V[num_rows];
   }
   return CalculateInnerProduct;
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
//    ScaleVector                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void ScaleVector(double *v, double x, int n)                  //
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
//     if ( x != 0.0)  ScaleVector(v, x,N);                       //
//      printf("The vector v is \n"); ...                                     //
////////////////////////////////////////////////////////////////////////////////
void ScaleVector(double X[], double TempScalar, unsigned int num_cols) 
{
   double z = 1.0 / TempScalar;
   for (; num_cols > 0; num_cols--) *X++ *= z;
}

std::vector<double> CreateVectorFilledWithValue(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
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

int Output_Data(vector<double> A[], double X[], double B[], unsigned int num_rows, unsigned int num_cols, double eigenvalue, double tolerance, unsigned int max_iter)
{
  unsigned int i;
  std::cout << std::setprecision (7) << endl;
  std::cout << "******************** Solve Ax = B ********************" << endl;
  // Displaying the 2D vector
  std::cout << "The vector A is the following: " << endl; 
  PrintMatrix(&A[0], num_cols);  
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

int main ()
{
   int max_iter {1000};
   double eigenvalue {0.0};
   double tolerance {0.00001};
   std::cout << std::setprecision(5);     
   int num_rows {0}, num_cols {0};
   matrix<double> A;
   ReadMatrixMarketMatrix(A,num_rows,num_cols);
   vector <double> B;
   ReadMatrixMarketVector(B, num_cols);      
   vector<double> X = CreateVectorFilledWithValue(num_rows);
   vector<double> array_A = CreateVectorFilledWithValue(num_rows*num_cols);
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
