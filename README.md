# Matrix Operations in C

## Overview

This C program provides a set of functions to perform various operations on matrices, including addition, subtraction, multiplication, inversion, and more. The program uses a custom-defined `Matrix` structure and supports operations on matrices of different orders.

## Features

- **Matrix Operations**: Addition, subtraction, multiplication, transpose, determinant, inversion, and more.
- **Symmetry Check**: Check if a matrix is symmetric.
- **Trace Calculation**: Calculate the trace of a square matrix.
- **Identity Matrix**: Generate an identity matrix of a given order.
- **Scalar Multiplication**: Multiply a matrix by a scalar.
- **Random Matrix Generation**: Generate a matrix with random values within a specified range.
- **Matrix Comparison**: Check if two matrices are equal.

## Getting Started

To use the matrix operations in your C program:

1. Include the `matrices.h` header file in your source file.
2. Create matrices using the `mat_alloc` function.
3. Perform various operations on matrices using the provided functions.

```c
#include "matrices.h"

int main() {
    // Example usage
    Matrix A = mat_alloc(3, 3);
    // ... perform matrix operations ...
    return 0;
}
```

## Functions

- `mat_alloc`: Allocate memory for a matrix.
- `mat_populate`: Populate a matrix with data.
- `isSquareMatrix`: Check if a matrix is square.
- `isEqual`: Check if two matrices are equal.
- `transpose`: Transpose a matrix.
- `isSym`: Check if a matrix is symmetric.
- `clearMatrix`: Set all elements of a matrix to 0.
- `viewMatrix`: Print a matrix.
- `randomizeMatrix`: Fill a matrix with random values.
- ... and more.


## Acknowledgments

- Author: Samarth Pyati
- Version: 1.1 (17/1/24)
