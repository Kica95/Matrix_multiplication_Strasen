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