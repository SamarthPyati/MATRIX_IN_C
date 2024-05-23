/**
    Matrix Multiplication
    matrices.c (improved branch)

    "Working with pass by references instead of pass by value
    to increase performance and lower memory usage."
    Matrix data structure in C.

    @author Samarth Pyati
    @version 1.3 12/2/24
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

/* Macros for Handling Error Messages - From man pages of pthread_*/
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
        fprintf(stderr, "MAT_ERR: %s\n", msg);      \
        exit(EXIT_FAILURE);                         \
    } while (0)


typedef struct
{
    /* THE MAIN STRUCTURE OF MATRIX */
    double *data;
    size_t rows;
    size_t cols;
} Matrix;


/* Fetch the element at a certain position (i, j) in the matrix*/
#define MAT_AT(m, i, j) (m).data[(i) * (m).cols + (j)]
#define MAT_AT_P(m, i, j) (m)->data[(i) * (m)->cols + (j)] // P -> pointer refr

/* Explicit declarations of positions of elements in the matrices*/
#define MAT_AT_EX(m, i, cols, j) (m).data[(i) * (cols) + (j)]
#define MAT_AT_EX_P(m, i, cols, j) (m)->data[(i) * (cols) + (j)]

void getMatrixOrder(Matrix *m)
{
    printf("%zu x %zu\n", m->rows, m->cols);
}

size_t getTotalElements(Matrix *m)
{
    return m->rows * m->cols;
}


Matrix *mat_alloc(size_t rows, size_t cols)
{
    /* Allocate certain amount of memory to matrix */
    Matrix *m = (Matrix *) malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    // Not casting types explicitly to work with diff types later
    m->data = malloc(sizeof(*m->data) * rows * cols);
    if (m->data == NULL)
        HANDLE_ERROR_MSG("MAT_ALLOC: Memory Allocation Failed");
    return m;
}

void free_mat(Matrix *m)
{
    /* Free the memory which has been mallocated to matrix */
    free(m->data);
    free(m);
}

void mat_populate(Matrix *m, double *data_)
{
    /* Populate a given array into the matrix */
    memcpy(m->data, data_, sizeof(double) * m->rows * m->cols);
}

int isSquareMatrix(Matrix *matrix)
{
    return (matrix->rows == matrix->cols);
}

int isEqual(Matrix *a, Matrix *b)
{
    /* Check equality of two given matrices */
    if (!(a->rows == b->rows && a->cols == b->cols))
    {
        return 0;
    }

    for (size_t i = 0; i < a->rows; ++i)
    {
        for (size_t j = 0; j < a->cols; ++j)
        {
            if (MAT_AT_P(a, i, j) != MAT_AT_P(b, i, j))
                return 0;
        }
    }
    return 1;
}

Matrix *transpose(Matrix *m);
int isSym(Matrix *m)
{
    /* Check if a matrix is symmentric */
    if (isEqual(m, transpose(m)))
        return 1;
    return 0;
}

void clear(Matrix *m)
{
    /* Set all elements in matrix to 0 */
    for (int i = 0; i < m->rows; ++i)
    {
        for (int j = 0; j < m->cols; ++j)
        {
            MAT_AT_P(m, i, j) = 0;
        }
    }
}

// row Major
void view(Matrix *m, unsigned viewOpt)
{
    /* Print the matrix in a formatted Matrix */
    printf("[");
    for (int i = 0; i < m->rows; i++)
    {
        printf("[");
        for (int j = 0; j < m->cols; j++)
        {
            switch (viewOpt)
            {
            case 1:
                printf(" %d ", (int)MAT_AT_P(m, i, j));
                break;
            case 2:
                printf(" %g ", MAT_AT_P(m, i, j));
                break;
            default:
                printf(" %f ", MAT_AT_P(m, i, j));
                break;
            }

            if (j < m->cols - 1)
                printf(","); // comma
        }
        printf("]");
        if (i < m->rows - 1)
            printf("\n"); // newline
    }
    printf("]\n\n");
}

void mat_rand(Matrix *m, const int MIN, const int MAX)
{
    /* Generate a random matrix with int elements in range(MIN, MAX) */
    if (MIN > MAX)
        HANDLE_ERROR_MSG("RAND: Incorrect assignment to Max and Min (MIN > MAX)");
    const int SIZE = m->rows * m->cols;
    for (unsigned int i = 0; i < SIZE; ++i)
    {
        m->data[i] = rand() % (MAX - MIN + 1) + MIN;
    }
}

void mat_randf(Matrix *m, const int MIN, const int MAX)
{
    /* Generate a random matrix with double elements in range(MIN, MAX) */
    if (MIN > MAX)
        HANDLE_ERROR_MSG("RAND_F: Incorrect assignment to Max and Min (MIN > MAX)");

    const int SIZE = m->rows * m->cols;

    for (unsigned int i = 0; i < SIZE; ++i)
    {
        double el = (double)rand() / (double)RAND_MAX; // float value between 0 and 1
        m->data[i] = el * (MAX - MIN) + MIN;
    }
}


Matrix *mat_add(Matrix *a, Matrix *b)
{
    /* Add two given matrices */
    if (!(a->rows == b->rows && a->cols == b->cols))
    {
        HANDLE_ERROR_MSG("ADD: Matrices should have same order in Addition");
    }

    Matrix *result = mat_alloc(a->rows, a->cols);

    for (size_t i = 0; i < a->rows; i++)
    {
        for (size_t j = 0; j < a->cols; j++)
        {
            MAT_AT_P(result, i, j) = MAT_AT_P(a, i, j) + MAT_AT_P(b, i, j);
        }
    }
    return result;
}

Matrix *mat_sub(Matrix *a, Matrix *b)
{
    /* Subtract two given matrices */
    if (!(a->rows == b->rows && a->cols == b->cols))
    {
        HANDLE_ERROR_MSG("SUB: Matrices should have same order for Subtraction");
    }

    Matrix *result = mat_alloc(a->rows, a->cols);

    for (size_t i = 0; i < a->rows; i++)
    {
        for (size_t j = 0; j < a->cols; j++)
        {
            MAT_AT_P(result, i, j) = MAT_AT_P(a, i, j) - MAT_AT_P(b, i, j);
        }
    }
    return result;
}

Matrix *mat_mul(Matrix *a, Matrix *b)
{
    /* Multiply Two given matrices */

    // Condition for mat_mul: A = 2 x 3 || B = 3 x 2, RESULT = 2 x 2

    if (a->cols != b->rows)
    {
        HANDLE_ERROR_MSG("MAT_MUL: Inappropriate order for matrix multiplication");
    }

    int rows = a->rows;
    int cols = a->cols;

    Matrix *result = mat_alloc(a->rows, b->cols);

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            // sum = a[i][k] + b[k][j]
            int sum = 0;
            for (size_t k = 0; k < cols; k++)
            {
                // sum += (*(a->data + i * cols + k)) * (*(b->data + k * cols + j));
                sum += MAT_AT_EX_P(a, i, cols, k) * MAT_AT_EX_P(b, k, cols, j);
            }
            // *(result->data + i * cols + j) = sum;
            MAT_AT_EX_P(result, i, result->cols, j) = sum;
        }
    }
    return result;
}

Matrix *transpose(Matrix *m)
{
    /* Calculate Transpose of Matrix */
    size_t rows = m->cols;
    size_t cols = m->rows;
    // A, A` | A[i][j] = A`[j][i]
    Matrix *result = mat_alloc(rows, cols);
    // row & cols are swapped

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            MAT_AT_P(result, i, j) = MAT_AT_P(m, j, i);
        }
    }
    return result;
}


Matrix *minor(Matrix *m, int r, int c)
{
    /* Gets the minor Matrix of a element in a matrix */
    if (r < 0 || r >= m->rows && c < 0 || c >= m->cols)
    {
        HANDLE_ERROR_MSG("MINOR: Given row or column are out of bounds of the dimensions of matrix.");
    }

    Matrix *minor_ = mat_alloc(m->rows - 1, m->cols - 1);

    int minorRows = 0, minorCols = 0;

    for (int i = 0; i < m->rows; i++)
    {
        if (i == r)
        {
            continue; // skipping the row of minor_ element
        }

        minorCols = 0; // Reset minor_ column for each row in the original m

        for (int j = 0; j < m->cols; j++)
        {
            if (j == c)
            {
                continue;
            }

            // copying the desired data
            // *(minor_->data + minorRows * minor_->cols + minorCols) = *(m->data + i * m->cols + j);
            MAT_AT_P(minor_, minorRows, minorCols) = MAT_AT_P(m, i, j);
            minorCols++;
        }
        minorRows++;
    }
    return minor_;
}

double det(Matrix *m)
{
    /* Calculate determinant of the given matrix */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("DET: Determinant only defined for square matrices (ORDER: n x n)");
    }

    // Calculate determinant for Matrix of Order = 2
    if (m->rows == 2 && m->cols == 2)
    {
        return MAT_AT_P(m, 0, 0) * MAT_AT_P(m, 1, 1) - MAT_AT_P(m, 0, 1) * MAT_AT_P(m, 1, 0);
    }

    double determinant = 0.0;
    Matrix *minorMatrix = mat_alloc(m->rows - 1, m->cols - 1);

    for (size_t i = 0; i < m->rows; ++i)
    {
        for (size_t minor_i = 0, k = 0; k < m->rows; ++k)
        {
            if (k != i)
            {
                for (size_t j = 1; j < m->cols; ++j)
                {
                    MAT_AT_P(minorMatrix, minor_i, j - 1) = MAT_AT_P(m, k, j);
                }
                ++minor_i;
            }
        }

        int sign = (i % 2 == 0) ? 1 : -1;
        determinant += sign * MAT_AT_P(m, i, 0) * det(minorMatrix);
    }

    // Free memory for the minorMatrix
    free_mat(minorMatrix);

    return determinant;
}

Matrix *mat_minor(Matrix *m)
{
    /* Returns a matrix containing minors of all elements of given Matrix */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("MINOR: Minor Operations valid on Square Matrices only");
    }

    Matrix *result = mat_alloc(m->rows, m->cols);

    for (size_t i = 0; i < m->rows; i++)
    {
        for (size_t j = 0; j < m->cols; j++)
        {
            Matrix *tempMinorMatrix = minor(m, i, j);
            double det_ = det(tempMinorMatrix);
            // *(result->data + i * result->cols + j) = det_;
            MAT_AT_P(result, i, j) = det_;
            free_mat(tempMinorMatrix); // Free memory allocated for the temporary minor matrix
        }
    }
    return result;
}

double cof(Matrix *m, size_t r, size_t c)
{
    /* Calculate cofactor of specific element a(r, c) in Matrix m */
    if (r < 0 || r >= m->rows || c < 0 || c >= m->cols)
    {
        HANDLE_ERROR_MSG("COF: Given row or column are out of bounds of the dimensions of Matrix. ");
    }

    double result;
    Matrix *part = minor(m, r, c);
    result = pow(-1, (r + c + 2)) * det(part); // + 2 is added as it is 0 based indexed
    free_mat(part);
    return result;
}

Matrix *mat_cof(Matrix *m)
{
    /* Returns a matrix containing cofactors of all elements of given Matrix */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("M_COF: Cofactors operations valid on Square Matrices only");
    }

    Matrix *cofactor_matrix = mat_alloc(m->rows, m->cols);

    if (cofactor_matrix->data == NULL)
    {
        HANDLE_ERROR_MSG("M_COF: Memory Allocation Failed for Cofactor Matrix");
    }

    for (size_t i = 0; i < m->rows; ++i)
    {
        for (size_t j = 0; j < m->cols; ++j)
        {
            double result = cof(m, i, j);
            // *(cofactor_matrix->data + i * m->cols + j) = result;
            MAT_AT_EX_P(cofactor_matrix, i, m->cols, j) = result;
        }
    }

    return cofactor_matrix;
}

Matrix *adj(Matrix *m)
{
    /* Return adjoint of the Matrix */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("ADJ: Adjoint operations valid on Square Matrices only");
    }

    Matrix *adjoint = mat_alloc(m->rows, m->cols);

    for (size_t i = 0; i < m->rows; i++)
    {
        for (size_t j = 0; j < m->cols; j++)
        {
            Matrix *minor_mat = minor(m, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            int det_ = det(minor_mat);
            // *(adjoint->data + j * adjoint->cols + i) = sign * det_;
            MAT_AT_P(adjoint, j, i) = sign * det_;
            free_mat(minor_mat);
        }
    }

    return adjoint;
}

Matrix *inv(Matrix *m)
{
    /* Return inverse of a matrix */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("INV: Inverse operation valid on Square Matrices only");
    }

    double det_ = det(m);

    if (det_ == 0.0)
    {
        HANDLE_ERROR_MSG("INV: Inverse does not exist for a singular matrix (determinant is zero)");
    }

    Matrix *adj_ = adj(m);
    Matrix *inverse_matrix = mat_alloc(m->rows, m->cols);

    for (size_t i = 0; i < m->rows; i++)
    {
        for (size_t j = 0; j < m->cols; j++)
        {
            // Avoid division by zero
            if (det_ != 0.0)
            {
                MAT_AT_P(inverse_matrix, i, j) = MAT_AT_P(adj_, i, j) / det_;
            }
            else
            {
                HANDLE_ERROR_MSG("INV: Division by zero in inverse calculation");
            }
        }
    }

    free_mat(adj_);

    return inverse_matrix;
}

Matrix *mabs(Matrix *m)
{
    /* Return absolute values of the elements */
    for (size_t i = 0; i < m->rows; ++i)
    {
        for (size_t j = 0; j < m->cols; ++j)
        {
            MAT_AT_P(m, i, j) = fabs(MAT_AT_P(m, i, j));
        }
    }
    return m;
}

double trace(Matrix *m)
{
    /* Return the trace of the matrix, i.e return sum of diagonal elements */
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("TRACE: Trace operations valid only on square Matrices.");
    }

    double res = 0;
    for (size_t i = 0; i < m->rows; ++i)
    {
        for (size_t j = 0; j < m->cols; ++j)
        {
            if (i == j)
                res += MAT_AT_P(m, i, j);
        }
    }
    return res;
}

Matrix *null_mat(size_t rows, size_t cols)
{
    /* Generate a null matrix */
    Matrix *null = mat_alloc(rows, cols);
    clear(null);
    return null;
}

Matrix *fill(size_t rows, size_t cols, double element)
{
    /* Generate a matrix with all entries as a specific element */
    Matrix *m = mat_alloc(rows, cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            MAT_AT_P(m, i, j) = element;
        }
    }
    return m;
}

Matrix *diag(double *arr, size_t SIZE)
{
    /* convert 1D arr to a diagonal matrix */
    Matrix *m = mat_alloc(SIZE, SIZE);
    for (size_t i = 0; i < m->rows; i++)
    {
        for (size_t j = 0; j < m->cols; j++)
        {
            if (i == j)
                MAT_AT_P(m, i, j) = arr[i];
        }
    }
    return m;
}

Matrix *identity(size_t SIZE)
{
    /* Generate an indentity matrix i.e, null matrix with diagonals as 1 */
    double arr[SIZE];
    for (unsigned int i = 0; i < SIZE; i++)
        arr[i] = 1;
    Matrix *id = diag(arr, SIZE);
    return id;
}

void mat_mul_k(Matrix *m, double k)
{
    /* Multiply a matrix with a scalar value */
    for (size_t i = 0; i < m->rows; i++)
    {
        for (size_t j = 0; j < m->cols; j++)
        {
            MAT_AT_P(m, i, j) *= k;
        }
    }
}

Matrix *sum(Matrix *m, unsigned int axis)
{
    /* Calculate sum of a matrix along a axis. */
    Matrix *result;
    switch (axis)
    {
        case 0: // whole sum
            result = mat_alloc(1, 1);
            double sum = 0.0;
            for (size_t i = 0; i < m->rows; i++)
            {
                for (size_t j = 0; j < m->cols; j++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
            }
            MAT_AT_P(result, 0, 0) = sum;
            break;
        case 1: // Row-wise
            result = mat_alloc(m->rows, 1);
            for (size_t i = 0; i < m->rows; i++)
            {
                double sum = 0.0;
                for (size_t j = 0; j < m->cols; j++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
                MAT_AT_P(result, i, 0) = sum;
            }
            break;
        case 2: // Column-wise
            result = mat_alloc(1, m->cols);
            for (size_t j = 0; j < m->cols; j++)
            {
                double sum = 0.0;
                for (size_t i = 0; i < m->rows; i++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
                MAT_AT_P(result, 0, j) = sum;
            }
            break;
        default:
            // Handle invalid axis value
            HANDLE_ERROR_MSG("SUM: Invalid axis value. Use 0 to get whole sum 1 for row-wise and 2 for column-wise.");
            break;
    }
    return result;
}

Matrix *mean(Matrix *m, unsigned int axis)
{
    /* Calculate mean of a matrix along a axis. */
    Matrix *result;

    switch (axis)
    {
        case 0: // whole mean
            result = mat_alloc(1, 1);
            double sum = 0.0;
            for (size_t i = 0; i < m->rows; i++)
            {
                for (size_t j = 0; j < m->cols; j++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
            }
            double mean = sum / getTotalElements(m);
            MAT_AT_P(result, 0, 0) = mean;
            break;
        case 2: // Row-wise
            result = mat_alloc(m->rows, 1);
            for (size_t i = 0; i < m->rows; i++)
            {
                double sum = 0.0;
                for (size_t j = 0; j < m->cols; j++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
                double mean = sum / m->cols;
                MAT_AT_P(result, i, 0) = mean;
            }
            break;
        case 1: // Column-wise
            result = mat_alloc(1, m->cols);
            for (size_t j = 0; j < m->cols; j++)
            {
                double sum = 0.0;
                for (size_t i = 0; i < m->rows; i++)
                {
                    sum += MAT_AT_P(m, i, j);
                }
                double mean = sum / m->rows;
                MAT_AT_P(result, 0, j) = mean;
            }
            break;
        default:
            // Handle invalid axis value
            HANDLE_ERROR_MSG("MEAN: Invalid axis value. Use 0 to get whole mean, 1 for row-wise mean, and 2 for column-wise mean.");
            break;
    }

    return result;
}

Matrix *dev(Matrix *m)
{
    /* Calculate Deviation of a matrix */
    Matrix *dev = mat_alloc(m->rows, m->cols);
    double _mean = MAT_AT_P(mean(m, 0), 0, 0);

    for (size_t i = 0; i < m->rows; ++i)
    {
        for (size_t j = 0; j < m->cols; ++j)
        {
            MAT_AT_P(dev, i, j) = MAT_AT_P(m, i, j) - _mean;
        }
    }
    return dev;
}


double mean_dev(Matrix *m)
{
    /* Calculate Mean Deviation of a matrix */
    Matrix *_devs = mabs(dev(m)); // absolute deviations
    return MAT_AT_P(mean(_devs, 0), 0, 0);
}


Matrix *var(Matrix *m, unsigned int axis)
{
    /* Calculate variance of a matrix along a axis. */
    Matrix *_dev = dev(m);
    Matrix *res;

    switch (axis)
    {
        case 0:
            res = mat_alloc(1, 1);
            double _res = 0;
            for (size_t i = 0; i < _dev->rows; ++i)
            {
                for (size_t j = 0; j < _dev->cols; ++j)
                {
                    _res += pow(MAT_AT_P(_dev, i, j), 2);
                }
            }
            _res /= (double)getTotalElements(_dev);
            MAT_AT_P(res, 0, 0) = _res;
            break;
        case 1:         // row wise
            res = mat_alloc(m->rows, 1);
            for (size_t i = 0; i < m->rows; ++i)
            {
                double _res = 0.0;
                for (size_t j = 0; j < m->cols; ++j)
                {
                    _res += pow(MAT_AT_P(_dev, i, j), 2);
                }
                _res /= (double) m->rows;
                MAT_AT_P(res, i, 0) = _res;
            }
            break;
        case 2:         // column wise
            res = mat_alloc(1, m->cols);
            for (size_t j = 0; j < m->cols; ++j)
            {
                double _res = 0.0;
                for (size_t i = 0; i < m->rows; ++i)
                {
                    _res += pow(MAT_AT_P(_dev, i, j), 2);
                }
                _res /= (double) m->cols;
                MAT_AT_P(res, 0, j) = _res;
            }
            break;
        default:
            HANDLE_ERROR_MSG("VAR: Invalid axis value. Use 0 to get whole variance, 1 for row-wise variance, and 2 for column-wise variance.");
            break;
    }

    return res;
}

Matrix *std(Matrix *m, unsigned int axis)
{
    /* Calculate standard deviation of a matrix along a axis. */
    Matrix *res;

    switch (axis)
    {
    case 0:
        res = mat_alloc(1, 1);
        MAT_AT_P(res, 0, 0) = sqrt(MAT_AT_P(var(m, 0), 0, 0));
        break;
    case 1: // row wise
        res = mat_alloc(m->rows, 1);
        for (size_t i = 0; i < m->rows; ++i)
        {
            MAT_AT_P(res, i, 0) = sqrt(MAT_AT_P(var(m, 1), i, 0));
        }
        break;
    case 2:
        res = mat_alloc(1, m->cols);
        for (size_t j = 0; j < m->cols; j++)
        {
            MAT_AT_P(res, 0, j) = sqrt(MAT_AT_P(var(m, 2), 0, j));
        }
        break;
    default:
        HANDLE_ERROR_MSG("STD: Invalid axis value. Use 0 to get whole standard deviation, 1 for row-wise and 2 for column-wise standard deviation.");
        break;
    }

    return res;
}

Matrix *aug(Matrix *a, Matrix *b)
{
    /* Generate an augmented matrix from two given matrices a & b -> [a | b] */
    if (a->rows != b->rows)
    {
        HANDLE_ERROR_MSG("AUG: Can`t augment matrices. Number of rows not same.");
    }

    Matrix *aug = mat_alloc(a->rows, (a->cols + b->cols));

    for (size_t i = 0; i < aug->rows; ++i)
    {
        for (size_t j = 0; j < aug->cols; ++j)
        {
            if (j < a->cols)
            {
                MAT_AT_P(aug, i, j) = MAT_AT_P(a, i, j);
            }
            else
            {
                MAT_AT_P(aug, i, j) = MAT_AT_P(b, i, j - a->cols);
            }
        }
    }
    return aug;
}


#endif // MATRIX_H
