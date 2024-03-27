#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <chrono>
using namespace std;
template<typename T>
using matrix = std::vector<std::vector<T>>;

int CalculateDominantEigenvalue (double array_A[], double B[], double X[], int num_rows,
		      int num_cols, double eigenvalue, double tolerance, int max_iter);
int Eigenvalue_Power_Method (double array_A[], int n, double *eigenvalue, double x[],
			     double x0[], double tolerance, int max_tries);
int WriteSolutionVector (double array[], double x0[], double x[], int num_rows, int num_cols, double tolerance, int max_iter);		
int Power_Eigenvalue_Input();
int ConvertMatrixToDenseArray(double array_A[], double X[], double B[], int num_rows, int num_cols);
double CalculateVectorNorm (double v[], int n);
double CalculateInnerProduct (double u[], double v[], int n);
void Multiply_Matrix_by_Vector (double u[], double *A, int num_rows, int num_cols, double v[]);
void ScaleVector (double v[], double x, int n);
int InitializeMatrixToZero (double *X, int num_rows, int num_cols);
void InitializeVectorToZero (double *X, int n);

int ReadMatrixMarketMatrix(matrix<double>& A, int& num_rows, int& num_cols);
int ReadMatrixMarketVector(vector<double>& B, int& num_cols);
int Write1DArrayToMatrixMarketFile(const double B[], int num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols);
int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols);
int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols);
void PrintVector(const std::vector<double>& vector);
void PrintMatrix(const matrix<double>& A); 
std::vector<double> CreateVectorFilledWithValue(int num_rows);

int CalculateDominantEigenvalue (double array_A[], double B[], double X[], int num_rows, int num_cols, double eigenvalue, double tolerance, int max_iter)
{
  fprintf (stdout, "\n");
   //by default all values are 0
  // Enter code to process data here.
  int tries = { 1 }, err = {0};
  printf ("Prog: Power_Solution.cpp\n\n\n");
  //Eigenvalue_Power_Method (&array_A[0], num_cols, &eigenvalue, &X[0], &B[0], tolerance, max_iter);
  if (tries >= max_iter)
    return -1;
  err = Eigenvalue_Power_Method (&array_A[0], num_cols, &eigenvalue, &X[0], &B[0], tolerance, max_iter);
  if (err == -1)
    printf (" The iteration failed to converge\n");
  else if (err == -2)
    printf (" The iteration vector vanished\n");
  else if (err == -3)
    printf (" The estimate of the dominant eigenvalue vanished\n");
  //WriteSolutionVector (&array_A[0], &B[0], &X[0], num_rows, num_cols, tolerance, max_iter);
  return 0;
}

int Eigenvalue_Power_Method (double array_A[], int n, double *eigenvalue, double x[],
			 double x0[], double tolerance, int max_iter)
{
  int tries;
  double old_estimate;
  double dominantEigenvalue;
  double *px = x, *px0 = x0, *pdum;
  // Check that the user specified tolerance is greater than the
  // machine epsilon.  If not, set it to the machine epsilon.
  if (tolerance < DBL_EPSILON)
    tolerance = DBL_EPSILON;
  // Check that the initial vector is non-zero.
  // eigenvalue = Vector_L2_Norm(px0, n);
  *eigenvalue = CalculateVectorNorm (px0, n);
  if (*eigenvalue == 0.0)
    return -2;
  // Normalize the initial vector.
  ScaleVector (px0, *eigenvalue, n);
  // Calculate initial estimate of the dominant eigenvalue.
  ConvertMatrixToDenseArray(&array_A[0], &px0[0], &px[0], n, n);
  *eigenvalue = CalculateInnerProduct (px, px0, n);
  // Iterate for a maximum of max_tries times using the power method.
  // Exit if an error occurs or if the relative difference of two
  // successive estimates is less than the specified tolerance.
  for (tries = 1; tries < max_iter; tries++)
    {
      old_estimate = *eigenvalue;
      pdum = px0;
      px0 = px;
      px = pdum;
      *eigenvalue = CalculateVectorNorm (px0, n);
      dominantEigenvalue = *eigenvalue;
      if (*eigenvalue == 0.0)
	return -2;
      ScaleVector (px0, *eigenvalue, n);
      ConvertMatrixToDenseArray(&array_A[0], &px0[0], &px[0], n, n);
      *eigenvalue = CalculateInnerProduct (px, px0, n);
      if (*eigenvalue == 0.0)
	return -3;
      if (fabs ((*eigenvalue - old_estimate) / *eigenvalue) < tolerance)
       {

         break;
       }       
    }
  if (tries >= max_iter)
  {
    return -1;
  }
  //printf("The eigenvalue is %lf.\n", *eigenvalue);
  printf("The dominantEigenvalue is %lf.\n", dominantEigenvalue);
  return 0;
}

int WriteSolutionVector (double array[], double x0[], double x[], int num_rows, int num_cols, double tolerance, int max_iter)
{
  int i, j;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("A =\n");
  for (i = 0; i < num_rows; i++)
    {
      for (j = 0; j < num_cols; j++)
	{
	  printf ("%6.4lf   ", array[i * num_rows + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("B =\n");
  for (i = 0; i < num_cols; i++)
    printf ("%6.4lf   ", x0[i]);
  printf ("\n\n");
  FILE *myFile3;
  myFile3 = fopen("X.dat","w+");
  if (myFile3 == NULL)
  {
    printf("Error writing to file.\n");
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_X.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, num_rows, num_cols);
  printf ("The solution vector X is = \n");
  for (i = 0; i < num_cols; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, x[i]);
    printf ("%6.4lf    ", x[i]);
  }
  fclose(myFile3);
  printf ("\n\n");
  printf ("\n\n\n");
  printf ("Tolerance: %6.4lf\n", tolerance);
  printf ("max_iter: %d\n", max_iter);
  return 0;
}

int Power_Eigenvalue_Input()
{	  
  int num_rows {0}, num_cols {0};
  double eigenvalue {0.0}, tolerance {0.00001};
  int max_iter {1000};
  std::cout << std::setprecision(5);
  matrix<double> A;
  ReadMatrixMarketMatrix(A,num_rows,num_cols);
  vector <double> B;
  ReadMatrixMarketVector(B, num_cols);      
  vector<double> X = CreateVectorFilledWithValue(num_rows);
  vector<double> array_A = CreateVectorFilledWithValue(num_rows*num_cols);
   // Copy A matrix into U and V arrays
   for (int i = 0; i < num_rows; ++i) {
     for (int j = 0; j < num_cols; ++j) {
       array_A[i * num_rows + j] = A[i][j];
     }
   }
   ReadMatrixMarketVector(B, num_cols);   
   // Using time point and system_clock
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();	 
   // Begin to place new code here.
   CalculateDominantEigenvalue (&array_A[0], &B[0], &X[0], num_rows, num_cols, eigenvalue, tolerance, max_iter);
   WriteSolutionVector (&array_A[0], &B[0], &X[0], num_rows, num_cols, tolerance, max_iter);  
   // End the new code here. 
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);
   std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
   // Memory allocation released here.         
   array_A.clear();
   array_A.shrink_to_fit();  	
   A.clear();
   A.shrink_to_fit();
   B.clear();
   B.shrink_to_fit();
   X.clear();
   X.shrink_to_fit();
  return 0;	
}

int ConvertMatrixToDenseArray(double array_A[], double X[], double B[], int num_rows, int num_cols)
{
   int i,j;
   double sum = 0.0;
   for (i = 0; i < num_rows; i++)
   {
     sum = 0.0;
     for (j = 0; j < num_cols; j++)
     {
       sum += array_A[i*num_cols+j] * X[j];
       B[i] = sum;
     }
   }
   return 0;
}

double CalculateVectorNorm (double v[], int n)
{
  double eigenvalue_temp = 0.0;
  double x;
  int i;
  for (i = 0; i < n; i++)
    if (eigenvalue_temp < (x = fabs (v[i])))
      eigenvalue_temp = x;
  return eigenvalue_temp;
}

double CalculateInnerProduct (double u[], double v[], int n)
{
  double CalculateInnerProduct = 0.0;
  for (n--; n >= 0; n--)
    CalculateInnerProduct += u[n] * v[n];
  return CalculateInnerProduct;
}
void Multiply_Matrix_by_Vector (double u[], double *A, int num_rows, int num_cols,
			   double v[])
{
  int i, j;
  for (i = 0; i < num_rows; A += num_cols, i++)
    for (u[i] = 0.0, j = 0; j < num_cols; j++)
      u[i] += A[j] * v[j];
}
void ScaleVector (double v[], double x, int n)
{
  double z = 1.0 / x;
  for (; n > 0; n--)
    *v++ *= z;
}

int InitializeMatrixToZero (double *X, int num_rows, int num_cols)
{
  int n = num_rows * num_cols;
  for (; n > 0; n--)
    *X++ = 0.0;
  return 0;
}
 
void InitializeVectorToZero (double *X, int n)
{
  for (; n > 0; n--)
    *X++ = 0.0;
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

std::vector<double> CreateVectorFilledWithValue(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
}

int main ()
{
  Power_Eigenvalue_Input();
   return 0;
}
