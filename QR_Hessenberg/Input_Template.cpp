#include "Input_Template.h"
int ReadTridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols)
{ 
  double *superdiagonal;
  superdiagonal = new double[num_cols];
  double *diagonal;
  diagonal = new double[num_cols];
  double *subdiagonal;
  subdiagonal = new double[num_cols];
  for (int i = 0; i < num_rows; i++)
  {
    for (int j = 0; j < num_cols; j++)
	{
	  if (i == j)
	  {
		diagonal[i] = A[i][j];
      }
      else if  (i == j-1)
      {
		superdiagonal[i] = A[i][j];
	  }
	  else if (i == j+1)
	  {
		subdiagonal[i] = A[i][j];
	  }
	}
  }
  //Tridiagonal (A, &subdiagonal[0], &diagonal[0], &superdiagonal[0], &B[0], X, num_rows, num_cols);
  return 0;
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

int WriteMatrixMarketFile(const matrix<double>& A, int num_rows, int num_cols) {
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

int WriteVectorMarketFile(const std::vector<double>& X, int num_rows, int num_cols) {
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

int WriteTridiagonalToVectorMarketFile(const matrix<double>& A, double B[], vector<double>& X, int num_rows, int num_cols)
{
  std::cout << endl << "Implementing the Tridiagonal method: " << std::endl;
  //printf ("A =\n");
  //PrintVector2D(A);  
  Write1DArrayToMatrixMarketFile(&B[0], num_rows);   
  //printf ("X =\n");
  //PrintVector(X);
  WriteVectorMarketFile(X, num_rows, num_cols);
  return 0;
}

void PrintVector(const std::vector<double>& vector) {
// std::cout << "Displaying the vector: " << endl;
 std::cout << std::setprecision(5);
 for (const double& value : vector) {
   std::cout << value << " ";
 }
 std::cout << std::endl;
}

void PrintVector2D(const matrix<double>& A) {
// std::cout << "Displaying the 2D vector:" << endl;
 for (const auto& row : A) {
   for (const double& value : row) {
     std::cout << value << " ";
   }
   std::cout << std::endl;
 }
}

std::vector<double> CreateZeroVector(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
}


