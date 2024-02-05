/**
    Matrix Multiplication
    matrices.c (improved Version)
    "Using 2D array instead of a singular Array."
    Matrix data structure in C.

    @author Samarth Pyati
    @version 1.2 28/1/24
*/

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#define True 1
#define False 0

/* Macros for Handling Error Messages */
#define HANDLE_ERROR(msg)   \
    do                      \
    {                       \
        perror(msg);        \
        exit(EXIT_FAILURE); \
    } while (0)

#define HANDLE_ERROR_EN(en, msg) \
    do                           \
    {                            \
        errno = en;              \
        perror(msg);             \
        exit(EXIT_FAILURE);      \
    } while (0)

#define HANDLE_ERROR_MSG(msg)                       \
    do                                              \
    {                                               \
        fprintf(stderr, "Matrix_Error: %s\n", msg); \
        exit(EXIT_FAILURE);                         \
    } while (0)

/* THE MAIN STRUCTURE OF MATRIX */
typedef struct
{
    double *data;
    size_t rows;
    size_t cols;
} Matrix;

/* Fetch the element at a certain position (i, j) in the matrix*/
#define MAT_AT(m, i, j) (m).data[(i) * (m).cols + (j)]
#define MAT_ATP(m, i, j) (m)->data[(i) * (m)->cols + (j)] // P -> pointer refr

/* explicit declarations of positions of elements in the matrices*/
#define MAT_AT_EX(m, i, cols, j) (m).data[(i) * (cols) + (j)]
#define MAT_AT_EX_P(m, i, cols, j) (m)->data[(i) * (cols) + (j)]

void getMatrixOrder(Matrix matrix)
{
    printf("%zu x %zu\n", matrix.rows, matrix.cols);
}

size_t getTotalElements(Matrix m)
{
    return m.rows * m.cols;
}

/* Allocate certain amount of memory to matrices */
Matrix mat_alloc(size_t rows, size_t cols)
{
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = malloc(sizeof(*m.data) * rows * cols);
    if (m.data == NULL)
        HANDLE_ERROR_MSG("Memory Allocation Failed");
    return m;
}

void mat_populate(Matrix *m, double *data_)
{
    memcpy(m->data, data_, sizeof(double) * m->rows * m->cols);
}

int isSquareMatrix(Matrix matrix)
{
    return (matrix.rows == matrix.cols);
}

int isEqual(Matrix a, Matrix b)
{
    // dimensions
    if (!(a.rows == b.rows && a.cols == b.cols))
    {
        return 0;
    }

    for (size_t i = 0; i < a.rows; ++i)
    {
        for (size_t j = 0; j < a.cols; ++j)
        {
            if (MAT_AT(a, i, j) != MAT_AT(b, i, j))
                return 0;
        }
    }
    return 1;
}

Matrix transpose(Matrix m);
int isSym(Matrix m)
{
    if (isEqual(m, transpose(m)))
        return 1;
    return 0;
}

void clear(Matrix matrix)
{
    for (int i = 0; i < matrix.rows; ++i)
    {
        for (int j = 0; j < matrix.cols; ++j)
        {
            MAT_AT(matrix, i, j) = 0;
        }
    }
}

// row Major
void view(Matrix matrix, unsigned viewOpt)
{
    printf("[");
    for (int i = 0; i < matrix.rows; i++)
    {
        printf("[");
        for (int j = 0; j < matrix.cols; j++)
        {
            switch (viewOpt)
            {
            case 1:
                printf(" %d ", (int)MAT_AT(matrix, i, j));
                break;
            case 2:
                printf(" %g ", MAT_AT(matrix, i, j));
                break;
            default:
                printf(" %f ", MAT_AT(matrix, i, j));
                break;
            }

            if (j < matrix.cols - 1)
                printf(","); // comma
        }
        printf("]");
        if (i < matrix.rows - 1)
            printf("\n"); // newline
    }
    printf("]\n\n");
}

void mat_rand(Matrix m, const int MIN, const int MAX)
{
    if (MIN > MAX)
        HANDLE_ERROR_MSG("Incorrect assignment to Max and Min (MIN > MAX)");
    const int SIZE = m.rows * m.cols;
    for (unsigned int i = 0; i < SIZE; ++i)
    {
        m.data[i] = rand() % (MAX - MIN + 1) + MIN;
    }
}

void mat_randf(Matrix m, const int MIN, const int MAX)
{
    if (MIN > MAX)
        HANDLE_ERROR_MSG("Incorrect assignment to Max and Min (MIN > MAX)");

    const int SIZE = m.rows * m.cols;

    for (unsigned int i = 0; i < SIZE; ++i)
    {
        double el = (double)rand() / (double)RAND_MAX; // float value between 0 and 1
        m.data[i] = el * (MAX - MIN) + MIN;
    }
}


Matrix mat_add(Matrix a, Matrix b)
{
    if (!(a.rows == b.rows && a.cols == b.cols))
    {
        HANDLE_ERROR_MSG("Matrices should have same order in Addition");
    }

    Matrix result = mat_alloc(a.rows, a.cols);

    for (size_t i = 0; i < a.rows; i++)
    {
        for (size_t j = 0; j < a.cols; j++)
        {
            // *(result.data + i * result.cols + j) = *(a.data + i * a.cols + j) + *(b.data + i * b.cols + j);
            MAT_AT(result, i, j) = MAT_AT(a, i, j) + MAT_AT(b, i, j);
        }
    }
    return result;
}

Matrix mat_sub(Matrix a, Matrix b)
{
    if (!(a.rows == b.rows && a.cols == b.cols))
    {
        HANDLE_ERROR_MSG("Matrices should have same order for Subtraction");
    }

    Matrix result = mat_alloc(a.rows, a.cols);

    for (size_t i = 0; i < a.rows; i++)
    {
        for (size_t j = 0; j < a.cols; j++)
        {
            // *(result.data + i * result.cols + j) = *(a.data + i * a.cols + j) - *(b.data + i * b.cols + j);
            MAT_AT(result, i, j) = MAT_AT(a, i, j) - MAT_AT(b, i, j);
        }
    }
    return result;
}

Matrix mat_mul(Matrix a, Matrix b)
{
    // A = 2 x 3 || B = 3 x 2, RESULT = 2 x 2
    if (a.cols != b.rows)
    {
        HANDLE_ERROR_MSG("Inappropriate order for matrix multiplication");
    }

    int rows = a.rows;
    int cols = a.cols;

    Matrix result = mat_alloc(a.rows, b.cols);

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            // sum = a[i][k] + b[k][j]
            int sum = 0;
            for (size_t k = 0; k < cols; k++)
            {
                // sum += (*(a.data + i * cols + k)) * (*(b.data + k * cols + j));
                sum += MAT_AT_EX(a, i, cols, k) * MAT_AT_EX(b, k, cols, j);
            }
            // *(result.data + i * cols + j) = sum;
            MAT_AT_EX(result, i, cols, j) = sum;
        }
    }
    return result;
}

Matrix transpose(Matrix m)
{
    int rows = m.cols;
    int cols = m.rows;
    // A, A` | A[i][j] = A`[j][i]
    Matrix result = mat_alloc(rows, cols);
    // row & cols are swapped

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            // *(result.data + i * cols + j) = *(m.data + j * rows + i);
            MAT_AT(result, i, j) = MAT_AT(m, j, i);
        }
    }
    return result;
}

/* Gets the minor Matrix of a element in a matrix */
Matrix minor(Matrix matrix, int r, int c)
{
    if (r < 0 || r >= matrix.rows && c < 0 || c >= matrix.cols)
    {
        HANDLE_ERROR_MSG("Given row or column are out of bounds of the dimensions of Matrix.");
    }

    Matrix minor_ = mat_alloc(matrix.rows - 1, matrix.cols - 1);

    int minorRows = 0, minorCols = 0;

    for (int i = 0; i < matrix.rows; i++)
    {
        if (i == r)
        {
            continue; // skipping the row of minor_ element
        }

        minorCols = 0; // Reset minor_ column for each row in the original matrix

        for (int j = 0; j < matrix.cols; j++)
        {
            if (j == c)
            {
                continue;
            }

            // copying the desired data
            // *(minor_.data + minorRows * minor_.cols + minorCols) = *(matrix.data + i * matrix.cols + j);
            MAT_AT(minor_, minorRows, minorCols) = MAT_AT(matrix, i, j);
            minorCols++;
        }
        minorRows++;
    }
    return minor_;
}

double det(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Determinant only defined for square matrices (ORDER: n x n)");
    }

    // Calculate determinant for Matrix of Order = 2
    if (m.rows == 2)
    {
        return MAT_AT(m, 0, 0) * MAT_AT(m, 1, 1) - MAT_AT(m, 0, 1) * MAT_AT(m, 1, 0);
    }

    double determinant = 0;
    Matrix minorMatrix = mat_alloc(m.rows - 1, m.cols - 1);

    for (size_t i = 0; i < m.rows; ++i)
    {
        for (size_t minor_i = 0, k = 0; k < m.rows; ++k)
        {
            if (k != i)
            {
                for (size_t j = 1; j < m.cols; ++j)
                {
                    MAT_AT(minorMatrix, minor_i, j - 1) = MAT_AT(m, k, j);
                }
                ++minor_i;
            }
        }

        int sign = (i % 2 == 0) ? 1 : -1;
        determinant += sign * MAT_AT(m, i, 0) * det(minorMatrix);
    }

    // Free memory for the minorMatrix
    free(minorMatrix.data);

    return determinant;
}

Matrix mat_minor(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Minor Operations valid on Square Matrices only");
    }

    Matrix result = mat_alloc(m.rows, m.cols);

    for (size_t i = 0; i < m.rows; i++)
    {
        for (size_t j = 0; j < m.cols; j++)
        {
            Matrix tempMinorMatrix = minor(m, i, j);
            double det_ = det(tempMinorMatrix);
            // *(result.data + i * result.cols + j) = det_;
            MAT_AT(result, i, j) = det_;
            free(tempMinorMatrix.data); // Free memory allocated for the temporary minor matrix
        }
    }
    return result;
}

double cof(Matrix m, size_t r, size_t c)
{
    if (r < 0 || r >= m.rows && c < 0 || c >= m.cols)
    {
        HANDLE_ERROR_MSG("Given row or column are out of bounds of the dimensions of Matrix. ");
    }

    double result;
    Matrix part = minor(m, r, c);
    result = pow(-1, (r + c + 2)) * det(part); // + 2 is added as it is 0 based indexed
    return result;
}

Matrix mat_cof(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Cofactors operations valid on Square Matrices only");
    }

    Matrix cofactor_matrix = {(double *)malloc(sizeof(double) * m.rows * m.cols), m.rows, m.cols};

    if (cofactor_matrix.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed for Cofactor Matrix");
    }

    for (size_t i = 0; i < m.rows; ++i)
    {
        for (size_t j = 0; j < m.cols; ++j)
        {
            double result = cof(m, i, j);
            // *(cofactor_matrix.data + i * m.cols + j) = result;
            MAT_AT_EX(cofactor_matrix, i, m.cols, j) = result;
        }
    }

    return cofactor_matrix;
}

Matrix adj(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Adjoint operations valid on Square Matrices only");
    }

    Matrix adjoint = mat_alloc(m.rows, m.cols);

    for (size_t i = 0; i < m.rows; i++)
    {
        for (size_t j = 0; j < m.cols; j++)
        {
            Matrix minor_mat = minor(m, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            int det_ = det(minor_mat);
            // *(adjoint.data + j * adjoint.cols + i) = sign * det_;
            MAT_AT(adjoint, j, i) = sign * det_;
            free(minor_mat.data);
        }
    }

    return adjoint;
}

Matrix inv(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Inverse operation valid on Square Matrices only");
    }

    double det_ = det(m);

    if (det_ == 0.0)
    {
        HANDLE_ERROR_MSG("Inverse does not exist for a singular matrix (determinant is zero)");
    }

    Matrix adj_ = adj(m);
    Matrix inverse_matrix = mat_alloc(m.rows, m.cols);

    for (size_t i = 0; i < m.rows; i++)
    {
        for (size_t j = 0; j < m.cols; j++)
        {
            // Avoid division by zero
            if (det_ != 0.0)
            {
                MAT_AT(inverse_matrix, i, j) = MAT_AT(adj_, i, j) / det_;
            }
            else
            {
                HANDLE_ERROR_MSG("Division by zero in inverse calculation");
            }
        }
    }

    free(adj_.data);

    return inverse_matrix;
}

Matrix mabs(Matrix m)
{
    /* return absolute values of the elements */
    for (size_t i = 0; i < m.rows; ++i)
    {
        for (size_t j = 0; j < m.cols; ++j)
        {
            MAT_AT(m, i, j) = fabs(MAT_AT(m, i, j));
        }
    }
    return m;
}

double trace(Matrix m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Trace operations valid only on square Matrices.");
    }

    double res = 0;
    for (size_t i = 0; i < m.rows; ++i)
    {
        for (size_t j = 0; j < m.cols; ++j)
        {
            if (i == j)
                res += MAT_AT(m, i, j);
        }
    }
    return res;
}

Matrix null_mat(size_t rows, size_t cols)
{
    Matrix null = mat_alloc(rows, cols);
    clear(null);
    return null;
}

Matrix fill(size_t rows, size_t cols, double element)
{
    Matrix m = mat_alloc(rows, cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            MAT_AT(m, i, j) = element;
        }
    }
    return m;
}

Matrix diag(double *arr, size_t SIZE)
{
    // convert 1D arr to a diagonal matrix
    Matrix m = mat_alloc(SIZE, SIZE);
    for (size_t i = 0; i < m.rows; i++)
    {
        for (size_t j = 0; j < m.cols; j++)
        {
            if (i == j)
                MAT_AT(m, i, j) = arr[i];
        }
    }
    return m;
}

Matrix identity(size_t SIZE)
{
    // Square matrix with diagonal elements as 1
    double arr[SIZE];
    for (unsigned int i = 0; i < SIZE; i++)
        arr[i] = 1;
    Matrix id = diag(arr, SIZE);
    return id;
}

void mat_mul_k(Matrix m, double k)
{
    // Multiply a matrix with a scalar value
    for (size_t i = 0; i < m.rows; i++)
    {
        for (size_t j = 0; j < m.cols; j++)
        {
            MAT_AT(m, i, j) *= k;
        }
    }
}

Matrix sum(Matrix m, unsigned int axis)
{
    Matrix result;
    switch (axis)
    {
    case 0: // whole sum
        result = mat_alloc(1, 1);
        double sum = 0.0;
        for (size_t i = 0; i < m.rows; i++)
        {
            for (size_t j = 0; j < m.cols; j++)
            {
                sum += MAT_AT(m, i, j);
            }
        }
        MAT_AT(result, 0, 0) = sum;
        break;
    case 1: // Row-wise
        result = mat_alloc(m.rows, 1);
        for (size_t i = 0; i < m.rows; i++)
        {
            double sum = 0.0;
            for (size_t j = 0; j < m.cols; j++)
            {
                sum += MAT_AT(m, i, j);
            }
            MAT_AT(result, i, 0) = sum;
        }
        break;
    case 2: // Column-wise
        result = mat_alloc(1, m.cols);
        for (size_t j = 0; j < m.cols; j++)
        {
            double sum = 0.0;
            for (size_t i = 0; i < m.rows; i++)
            {
                sum += MAT_AT(m, i, j);
            }
            MAT_AT(result, 0, j) = sum;
        }
        break;
    default:
        // Handle invalid axis value
        HANDLE_ERROR_MSG("Invalid axis value. Use 0 to get whole sum 1 for row-wise and 2 for column-wise.");
        break;
    }
    return result;
}

Matrix mean(Matrix m, unsigned int axis)
{
    Matrix result;

    switch (axis)
    {
    case 0: // whole mean
        result = mat_alloc(1, 1);
        double sum = 0.0;
        for (size_t i = 0; i < m.rows; i++)
        {
            for (size_t j = 0; j < m.cols; j++)
            {
                sum += MAT_AT(m, i, j);
            }
        }
        double mean = sum / getTotalElements(m);
        MAT_AT(result, 0, 0) = mean;
        break;
    case 2: // Row-wise
        result = mat_alloc(m.rows, 1);
        for (size_t i = 0; i < m.rows; i++)
        {
            double sum = 0.0;
            for (size_t j = 0; j < m.cols; j++)
            {
                sum += MAT_AT(m, i, j);
            }
            double mean = sum / m.cols;
            MAT_AT(result, i, 0) = mean;
        }
        break;
    case 1: // Column-wise
        result = mat_alloc(1, m.cols);
        for (size_t j = 0; j < m.cols; j++)
        {
            double sum = 0.0;
            for (size_t i = 0; i < m.rows; i++)
            {
                sum += MAT_AT(m, i, j);
            }
            double mean = sum / m.rows;
            MAT_AT(result, 0, j) = mean;
        }
        break;
    default:
        // Handle invalid axis value
        HANDLE_ERROR_MSG("Invalid axis value. Use 0 to get whole mean, 1 for row-wise mean, and 2 for column-wise mean.");
        break;
    }

    return result;
}

Matrix dev(Matrix m)
{
    Matrix dev = mat_alloc(m.rows, m.cols);
    double _mean = MAT_AT(mean(m, 0), 0, 0);

    for (size_t i = 0; i < m.rows; ++i)
    {
        for (size_t j = 0; j < m.cols; ++j)
        {
            MAT_AT(dev, i, j) = MAT_AT(m, i, j) - _mean;
        }
    }
    return dev;
}


double mean_dev(Matrix m)
{
    Matrix _devs = mabs(dev(m)); // absolute deviations
    return MAT_AT(mean(_devs, 0), 0, 0);
}


Matrix var(Matrix m, unsigned int axis)
{
    Matrix _dev = dev(m);
    Matrix res;

    switch (axis)
    {
    case 0:
        res = mat_alloc(1, 1);
        double _res = 0;
        for (size_t i = 0; i < _dev.rows; ++i)
        {
            for (size_t j = 0; j < _dev.cols; ++j)
            {
                _res += pow(MAT_AT(_dev, i, j), 2);
            }
        }
        _res /= (double)getTotalElements(_dev);
        MAT_AT(res, 0, 0) = _res;
        break;
    case 1:         // row wise
        res = mat_alloc(m.rows, 1);
        for (size_t i = 0; i < m.rows; ++i)
        {
            double _res = 0.0;
            for (size_t j = 0; j < m.cols; ++j)
            {
                _res += pow(MAT_AT(_dev, i, j), 2);
            }
            _res /= (double) m.rows;
            MAT_AT(res, i, 0) = _res;
        }
        break;
    case 2:         // column wise
        res = mat_alloc(1, m.cols);
        for (size_t j = 0; j < m.cols; ++j)
        {
            double _res = 0.0;
            for (size_t i = 0; i < m.rows; ++i)
            {
                _res += pow(MAT_AT(_dev, i, j), 2);
            }
            _res /= (double) m.cols;
            MAT_AT(res, 0, j) = _res;
        }
        break;
    default:
        HANDLE_ERROR_MSG("Invalid axis value. Use 0 to get whole variance, 1 for row-wise variance, and 2 for column-wise variance.");
        break;
    }

    return res;
}

Matrix std(Matrix m, unsigned int axis)
{
    Matrix res;

    switch (axis)
    {
    case 0:
        res = mat_alloc(1, 1);
        MAT_AT(res, 0, 0) = sqrt(MAT_AT(var(m, 0), 0, 0));
        break;
    case 1: // row wise
        res = mat_alloc(m.rows, 1);
        for (size_t i = 0; i < m.rows; ++i)
        {
            MAT_AT(res, i, 0) = sqrt(MAT_AT(var(m, 1), i, 0));
        }
        break;
    case 2:
        res = mat_alloc(1, m.cols);
        for (size_t j = 0; j < m.cols; j++)
        {
            MAT_AT(res, 0, j) = sqrt(MAT_AT(var(m, 2), 0, j));
        }
        break;
    default:
        HANDLE_ERROR_MSG("Invalid axis value. Use 0 to get whole standard deviation, 1 for row-wise and 2 for column-wise standard deviation.");
        break;
    }

    return res;
}


#endif // MATRIX_H

// #ifdef MATRIX_IMPLEMENTATION
// // CODE
// #endif // MATRIX_IMPLEMENTATION
