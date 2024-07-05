// StrassensMethodForMat.Mul.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <intrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int isPowerOfTwo(int n)
{
	// Check if a number is a power of two
	return (n > 0) && ((n & (n - 1)) == 0);
}
//The expression(n & (n - 1)) removes the least significant set bit of n.If the result is 0, it means that n had only one bit set, indicating that it is a power of two.
//The condition(n > 0) ensures that we are only considering positive integers, as 0 is not a power of 2.
///////////////////////////////////////////////////////////////////////////////////////////////////////////
double** createRandomMatrix(int rows, int cols)
{
	double** matrix = (double**)malloc(rows * sizeof(double*));

	for (int i = 0; i < rows; i++)
	{
		matrix[i] = (double*)malloc(cols * sizeof(double));
		for (int j = 0; j < cols; j++)
		{
			matrix[i][j] = ((double)rand() / RAND_MAX) * 10.0; // Random double between 0 and 10
		}
	}

	return matrix;
}

// Function to print a matrix
void printMatrix(double** matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			printf("%8.3f ", matrix[i][j]);
		}
		printf("\n");
	}
}

// Function to transpose a matrix
double** transposeMatrix(double** matrix, int rows, int cols)
{
	double** transpose = (double**)malloc(cols * sizeof(double*));

	for (int i = 0; i < cols; i++)
	{
		transpose[i] = (double*)malloc(rows * sizeof(double));
		for (int j = 0; j < rows; j++)
		{
			transpose[i][j] = matrix[j][i];
		}
	}

	return transpose;
}

// Function to free memory allocated for a matrix
void freeMatrix(double** matrix, int rows)
{
	for (int i = 0; i < rows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

// Function to check and pad a matrix to make it square
double** makeSquareMatrix(double** matrix, int rows, int cols)
{
	// Determine the size of the new square matrix
	int NewSize = (rows > cols) ? rows : cols;// findin max dimension

	// Check if n is already a power of two
	if (isPowerOfTwo(NewSize))
	{
		printf("%d can be represented as 2^i.\n", NewSize);
	}
	else
	{
		// Find the next power of two
		while (!isPowerOfTwo(NewSize))
		{
			NewSize++;
		}
	}


	// Allocate memory for the new square matrix
	double** squareMatrix = (double**)malloc(NewSize * sizeof(double*));
	for (int i = 0; i < NewSize; ++i)
	{
		squareMatrix[i] = (double*)calloc(NewSize, sizeof(double));
	}

	// Copy the original matrix to the top-left part of the new square matrix
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			squareMatrix[i][j] = matrix[i][j];
		}
	}

	return squareMatrix;
}

// function to sum two matrices
void sum(double** a, double** b, double** result, int tam)
{

	int i, j;

	for (i = 0; i < tam; i++)
	{
		for (j = 0; j < tam; j++)
		{
			result[i][j] = a[i][j] + b[i][j];
		}
	}
}

/*------------------------------------------------------------------------------*/
// function to subtract two matrices
void subtract(double** a, double** b, double** result, int tam)
{

	int i, j;

	for (i = 0; i < tam; i++)
	{
		for (j = 0; j < tam; j++)
		{
			result[i][j] = a[i][j] - b[i][j];
		}
	}
}

/*------------------------------------------------------------------------------*/
// This function allocates the matrix using malloc, and initializes it. If the variable random is passed
// as zero, it initializes the matrix with zero, if it's passed as 1, it initializes the matrix with random
// values. If it is passed with any other int value (like -1 for example) the matrix is initialized with no
// values in it. The variable tam defines the length of the matrix.
double** allocate_real_matrix(int tam)
{

	int i, j, n = tam, m = tam;
	//	double** v;         // pointer to the vector

	// allocates one vector of vectors (matrix)
	double** v = (double**)malloc(n * sizeof(double*));

	if (v == NULL)
	{
		printf("** Error in matrix allocation: insufficient memory **");
		return (NULL);
	}

	// allocates each row of the matrix
	for (i = 0; i < n; i++)
	{
		v[i] = (double*)malloc(m * sizeof(double));

		if (v[i] == NULL)
		{
			printf("** Error: Insufficient memory **");
			freeMatrix(v, n);
			return (NULL);
		}

		for (j = 0; j < m; j++)
			v[i][j] = 0.0;//initalize to zero

	}

	return (v);     // returns the pointer to the vector. 
}
//////////////////////Strassen Method Function\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

 void strassen(double** a, double** b, double** c, int tam)
{


	if (tam == 2)
	{
		c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
		c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
		c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
		c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];
		return;
	}

	// other cases are treated here:
	else
	{
		int newTam = tam / 2;
		double** a11, ** a12, ** a21, ** a22;
		double** b11, ** b12, ** b21, ** b22;
		double** c11, ** c12, ** c21, ** c22;
		double** p1, ** p2, ** p3, ** p4, ** p5, ** p6, ** p7;

		// memory allocation:
		a11 = allocate_real_matrix(newTam);
		a12 = allocate_real_matrix(newTam);
		a21 = allocate_real_matrix(newTam);
		a22 = allocate_real_matrix(newTam);

		b11 = allocate_real_matrix(newTam);
		b12 = allocate_real_matrix(newTam);
		b21 = allocate_real_matrix(newTam);
		b22 = allocate_real_matrix(newTam);

		c11 = allocate_real_matrix(newTam);
		c12 = allocate_real_matrix(newTam);
		c21 = allocate_real_matrix(newTam);
		c22 = allocate_real_matrix(newTam);

		p1 = allocate_real_matrix(newTam);
		p2 = allocate_real_matrix(newTam);
		p3 = allocate_real_matrix(newTam);
		p4 = allocate_real_matrix(newTam);
		p5 = allocate_real_matrix(newTam);
		p6 = allocate_real_matrix(newTam);
		p7 = allocate_real_matrix(newTam);

		double** aResult = allocate_real_matrix(newTam);
		double** bResult = allocate_real_matrix(newTam);

		int i, j;

		//dividing the matrices in 4 sub-matrices:
		for (i = 0; i < newTam; i++)
		{
			for (j = 0; j < newTam; j++)
			{
				a11[i][j] = a[i][j];
				a12[i][j] = a[i][j + newTam];
				a21[i][j] = a[i + newTam][j];
				a22[i][j] = a[i + newTam][j + newTam];

				b11[i][j] = b[i][j];
				b12[i][j] = b[i][j + newTam];
				b21[i][j] = b[i + newTam][j];
				b22[i][j] = b[i + newTam][j + newTam];
			}
		}

		// Calculating p1 to p7:

		sum(a11, a22, aResult, newTam); // a11 + a22
		sum(b11, b22, bResult, newTam); // b11 + b22
		strassen(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)

		sum(a21, a22, aResult, newTam); // a21 + a22
		strassen(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)

		subtract(b12, b22, bResult, newTam); // b12 - b22
		strassen(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)

		subtract(b21, b11, bResult, newTam); // b21 - b11
		strassen(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)

		sum(a11, a12, aResult, newTam); // a11 + a12
		strassen(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)	

		subtract(a21, a11, aResult, newTam); // a21 - a11
		sum(b11, b12, bResult, newTam); // b11 + b12
		strassen(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)

		subtract(a12, a22, aResult, newTam); // a12 - a22
		sum(b21, b22, bResult, newTam); // b21 + b22
		strassen(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)

		// calculating c21, c21, c11 e c22:

		sum(p3, p5, c12, newTam); // c12 = p3 + p5
		sum(p2, p4, c21, newTam); // c21 = p2 + p4

		sum(p1, p4, aResult, newTam); // p1 + p4
		sum(aResult, p7, bResult, newTam); // p1 + p4 + p7
		subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7

		sum(p1, p3, aResult, newTam); // p1 + p3
		sum(aResult, p6, bResult, newTam); // p1 + p3 + p6
		subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6

		// Grouping the results obtained in a single matrix:
		for (i = 0; i < newTam; i++)
		{
			for (j = 0; j < newTam; j++)
			{
				c[i][j] = c11[i][j];
				c[i][j + newTam] = c12[i][j];
				c[i + newTam][j] = c21[i][j];
				c[i + newTam][j + newTam] = c22[i][j];
			}
		}

		// deallocating memory (free):
		freeMatrix(a11, newTam);
		freeMatrix(a12, newTam);
		freeMatrix(a21, newTam);
		freeMatrix(a22, newTam);

		freeMatrix(b11, newTam);
		freeMatrix(b12, newTam);
		freeMatrix(b21, newTam);
		freeMatrix(b22, newTam);

		freeMatrix(c11, newTam);
		freeMatrix(c12, newTam);
		freeMatrix(c21, newTam);
		freeMatrix(c22, newTam);

		freeMatrix(p1, newTam);
		freeMatrix(p2, newTam);
		freeMatrix(p3, newTam);
		freeMatrix(p4, newTam);
		freeMatrix(p5, newTam);
		freeMatrix(p6, newTam);
		freeMatrix(p7, newTam);

		freeMatrix(aResult, newTam);
		freeMatrix(bResult, newTam);
	} // end of else

} // end of Strassen function
////////////////////////////////Function for bit comparison for calculating nearest power of two dimenson  /////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{

	int m = 3;//rows
	int n = 3;//colunms

	int Dimension = (m > n) ? m : n;// findin max dimension

	// Check if n is already a power of two
	if (isPowerOfTwo(Dimension))
	{
		printf("%d can be represented as 2^i.\n", Dimension);
	}
	else
	{
		// Find the next power of two
		while (!isPowerOfTwo(Dimension))
		{
			Dimension++;
		}
	}


	double** matrix = createRandomMatrix(n, m);



	double** transposedMatrix = transposeMatrix(matrix, n, m);



	double** Rezult = allocate_real_matrix(Dimension);


	clock_t beginsim = clock();

	double** matrix_corrected = makeSquareMatrix(matrix, n, m);


	double** transposedMatrixCorrected = makeSquareMatrix(transposedMatrix, m, n);


	strassen(matrix_corrected, transposedMatrixCorrected, Rezult, Dimension);

	clock_t endsim = clock();

	printf("%Lf", (long double)(endsim - beginsim));


#pragma region Check
	printf("\nOriginal Matrix:\n");
	printMatrix(matrix, n, m);

	printf("\n Matrix Ttransposed:\n");
	printMatrix(transposedMatrix, m, n);

	printf("\n Matrix square:\n");
	printMatrix(matrix_corrected, Dimension, Dimension);

	printf("\nTransposed Matrix square:\n");
	printMatrix(transposedMatrixCorrected, Dimension, Dimension);

	printf("\n  Strassen method rezult :\n");
	printMatrix(Rezult, Dimension, Dimension);

	printf("\n rezult:\n");
	printMatrix(Rezult, n, n);

#pragma endregion


	freeMatrix(matrix, n);
	freeMatrix(transposedMatrix, m);
	freeMatrix(matrix_corrected, Dimension);
	freeMatrix(transposedMatrixCorrected, Dimension);
	freeMatrix(Rezult, Dimension);


	return 0;
}


