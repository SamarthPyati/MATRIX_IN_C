# MATRIX_IN_C

## Overview

This C program provides a comprehensive set of functions to perform various operations on matrices. The program uses a custom-defined `Matrix` structure and supports operations on matrices of different orders.

## Features

- **Matrix Operations**: Addition, subtraction, multiplication, transpose, determinant, inversion, and more.
- **Symmetry Check**: Determine if a matrix is symmetric.
- **Trace Calculation**: Calculate the trace of a square matrix.
- **Identity Matrix**: Generate an identity matrix of a given order.
- **Scalar Multiplication**: Multiply a matrix by a scalar.
- **Random Matrix Generation**: Generate a matrix with random values within a specified range.
- **Minor, Cofactor, and Adjoint Matrices**: Functions for minor, cofactor, and adjoint operations.

## Getting Started

To use the matrix operations in your C program:

1. Include the `matrix.h` header file in your source file.
2. Create matrices using the `mat_alloc` function.
3. Set the random seed generator to the current time.
4. View the matrix using `view` function.
5. Perform various operations on matrices using the provided functions.

```c
#include "matrix.h"

int main(void)
{
    // Example usage
    time_t t;               // seed Randomise function
    srand(t);

    // Allocate matrices and randomise the matrix
    Matrix *c = mat_alloc(3, 3);
    mat_randf(c, 1, 10);
    view(c, 2);

    Matrix *d = mat_alloc(3, 3);
    mat_randf(d, 10, 100);
    view(d, 2);

    // Operations on matrices
    view(aug(c, d), 2);
    view(identity(4), 1);
    view(var(c, 2), 0);
    view(mat_mul(c, d), 1);
    view(mat_add(c, d), 1);
    view(mat_sub(c, d), 1);
    view(mat_cof(c), 1);
    view(mat_mul(inv(d), d), 0);
    view(sum(c, 2), 1);
    view(mean(c, 1), 1);
    view(std(c, 2), 0);
    view(var(c, 0), 0);

    printf("Mean Deviation: %lf\n", mean_dev(c));
    printf("Det(A) & Det(B): %f , %f \n", det(c), det(d));
    printf("SQUARE ? %d\n", isSquareMatrix(c));
    printf("SYM ? %d\n", isSym(c));
    return 0;
}
```

## Functions

| Function Name | Description |
| ------------- | ----------- |
| `mat_alloc` | Allocate memory for a matrix. |
| `mat_populate` | Populate a matrix with data. |
| `mat_rand` | Randomizes values in a matrix within a range. |
| `mat_add`, `mat_sub`, `mat_mul` | Arithmetic operations for matrices. |
| `isSquareMatrix` | Check if a matrix is square. |
| `isEqual` | Check if two matrices are equal. |
| `transpose` | Transpose a matrix. |
| `isSym` | Check if a matrix is symmetric. |
| `clear` | Set all elements of a matrix to 0. |
| `view` | Print a matrix. |
| `randomize` | Fill a matrix with random values. |
| `det` | Calculate the determinant of a square matrix. |
| `mat_minor` | Generate a matrix of minors for each element. |
| `cof` and `mat_cof` | Calculate the cofactor of a single element and the cofactor matrix. |
| `adj` | Calculate the adjoint of a matrix. |
| `inv` | Calculate the inverse of a square matrix. |
| `trace` | Calculate the trace of a square matrix. |
| `null_mat` | Generate a matrix filled with zeros. |
| `fill` | Fill a matrix with a specified element. |
| `diag` | Create a diagonal matrix from a given array. |
| `identity` | Generate an identity matrix of a given order. |
| `mat_mul_k` | Multiply a matrix by a scalar. |
| `sum` | Calculate the sum of all elements in a matrix. |
| `mean` | Calculate mean across different axes. |
| `mean_dev` | Calculate mean deviation of a matrix. |

## Error Handling

Macros (`HANDLE_ERROR`, `HANDLE_ERROR_EN`, and `HANDLE_ERROR_MSG`) efficiently handle errors during matrix operations.

## Acknowledgments

- Author: Samarth Sanjay Pyati
- Version: 1.3 (12/2/24)

