/**
    Matrix Multiplication
    matrices.c (improved Version)
    "Using 2D array instead of a singular Array."
    Matrix data structure in C.

    @author Samarth Pyati
    @version 1.1 17/1/24
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
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

typedef struct
{
    double *data;
    int rows;
    int cols;
} Matrix;

void getMatrixOrder(Matrix *matrix)
{
    printf("%d x %d\n", matrix->rows, matrix->cols);
}

unsigned int getTotalElements(Matrix *m)
{
    return m->rows * m->cols;
}

int getMatrixElement(Matrix *m, unsigned int row, unsigned int col)
{
    return *(m->data + row * m->cols + col);
}

int isSquareMatrix(Matrix *matrix)
{
    return (matrix->rows == matrix->cols);
}

void clearMatrix(Matrix *matrix)
{
    for (int i = 0; i < matrix->rows; ++i)
    {
        for (int j = 0; j < matrix->cols; ++j)
        {
            *(matrix->data + i * matrix->cols + j) = 0;
        }
    }
}

// row Major
void viewMatrix(Matrix *matrix, int addSpace)
{
    printf("[");
    for (int i = 0; i < matrix->rows; i++)
    {
        printf("[");
        for (int j = 0; j < matrix->cols; j++)
        {
            double el = *(matrix->data + i * matrix->cols + j);
            addSpace ? printf(" %f ", el) : printf("%f", el);

            if (j < matrix->cols - 1)
                printf(","); // comma
        }
        printf("]");
        if (i < matrix->rows - 1)
            printf("\n"); // newline
    }
    printf("]\n\n");
}

void randomizeMatrix(Matrix *matrix, const int MAX_RANGE, const int MIN_RANGE)
{
    const int SIZE = matrix->rows * matrix->cols;
    for (unsigned int i = 0; i < SIZE; ++i)
    {
        matrix->data[i] = rand() % (MAX_RANGE - MIN_RANGE + 1) + MIN_RANGE;
    }
}

Matrix addMatrix(Matrix *a, Matrix *b)
{
    if (!(a->rows == b->rows && a->cols == b->cols))
    {
        HANDLE_ERROR_MSG("Matrices should have same order in Addition");
    }

    Matrix result = {(double *)malloc(sizeof(double) * a->rows * a->cols), a->rows, a->cols};

    if (result.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed");
    }

    for (unsigned int i = 0; i < a->rows; i++)
    {
        for (unsigned int j = 0; j < a->cols; j++)
        {
            *(result.data + i * result.cols + j) = *(a->data + i * a->cols + j) + *(b->data + i * b->cols + j);
        }
    }
    return result;
}

Matrix subtractMatrix(Matrix *a, Matrix *b)
{
    if (!(a->rows == b->rows && a->cols == b->cols))
    {
        HANDLE_ERROR_MSG("Matrices should have same order for Subtraction");
    }

    Matrix result = {(double *)malloc(sizeof(double) * a->rows * a->cols), a->rows, a->cols};

    if (result.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed");
    }

    for (unsigned int i = 0; i < a->rows; i++)
    {
        for (unsigned int j = 0; j < a->cols; j++)
        {
            *(result.data + i * result.cols + j) = *(a->data + i * a->cols + j) - *(b->data + i * b->cols + j);
        }
    }
    return result;
}

Matrix multiplyMatrix(Matrix *a, Matrix *b)
{
    // A = 2 x 3 || B = 3 x 2, ORDER = 2 x 2
    if (a->rows != b->cols)
    {
        HANDLE_ERROR_MSG("Inappropriate order for matrix multiplication");
    }

    int rows = a->rows;
    int cols = a->cols;
    Matrix result = {malloc(sizeof(double) * a->rows * b->cols), a->rows, b->cols};

    for (unsigned int i = 0; i < rows; i++)
    {
        for (unsigned int j = 0; j < cols; j++)
        {
            // sum = a[i][k] + b[k][j]
            int sum = 0;
            for (unsigned int k = 0; k < cols; k++)
            {
                sum += (*(a->data + i * cols + k)) * (*(b->data + k * cols + j));
            }
            *(result.data + i * cols + j) = sum;
        }
    }
    return result;
}

Matrix transpose(Matrix *m)
{
    int rows = m->cols;
    int cols = m->rows;
    // A, A` | A[i][j] = A`[j][i]
    Matrix result = {malloc(sizeof(double) * rows * cols), rows, cols};
    // row & cols are swapped

    for (unsigned int i = 0; i < rows; ++i)
    {
        for (unsigned int j = 0; j < cols; ++j)
        {
            *(result.data + i * cols + j) = *(m->data + j * rows + i);
        }
    }
    return result;
}

/* Gets the minor Matrix of a element in a matrix */
Matrix minor(Matrix *matrix, int r, int c)
{
    if (r < 0 || r >= matrix->rows && c < 0 || c >= matrix->cols)
    {
        HANDLE_ERROR_MSG("Given row or column are out of bounds of the dimensions of Matrix.");
    }

    Matrix minor_ = {malloc(sizeof(double) * (matrix->rows - 1) * (matrix->cols - 1)), matrix->rows - 1, matrix->cols - 1};

    if (minor_.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed For Matrix");
    }

    int minorRows = 0, minorCols = 0;

    for (int i = 0; i < matrix->rows; i++)
    {
        if (i == r)
        {
            continue; // skipping the row of minor_ element
        }

        minorCols = 0; // Reset minor_ column for each row in the original matrix

        for (int j = 0; j < matrix->cols; j++)
        {
            if (j == c)
            {
                continue;
            }

            // copying the desired data
            *(minor_.data + minorRows * minor_.cols + minorCols) = *(matrix->data + i * matrix->cols + j);
            minorCols++;
        }
        minorRows++;
    }
    return minor_;
}

double determinant(Matrix *m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Determinant only defined for square matrices (ORDER: n x n)");
    }

    // Calculate determinant for Matrix of Order = 2
    if (m->rows == 2)
    {
        return m->data[0] * m->data[3] - m->data[1] * m->data[2];
    }

    double det = 0;

    for (unsigned int i = 0; i < m->rows; ++i)
    {
        Matrix minorMatrix = minor(m, 0, i);
        int sign = (i % 2 == 0) ? 1 : -1;
        det += sign * m->data[i] * determinant(&minorMatrix);
        free(minorMatrix.data);
    }
    return det;
}

Matrix minorMatrix(Matrix *m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Minor Operations valid on Square Matrices only");
    }

    Matrix result = {malloc(sizeof(double) * m->rows * m->cols), m->rows, m->cols};

    if (result.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed for Matrix.");
    }

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            Matrix tempMinorMatrix = minor(m, i, j);
            double det = determinant(&tempMinorMatrix);
            *(result.data + i * result.cols + j) = det;
            free(tempMinorMatrix.data); // Free memory allocated for the temporary minor matrix
        }
    }
    return result;
}

double cofactor(Matrix *m, unsigned int r, unsigned int c)
{
    if (r < 0 || r >= m->rows && c < 0 || c >= m->cols)
    {
        HANDLE_ERROR_MSG("Given row or column are out of bounds of the dimensions of Matrix. ");
    }

    double result;
    Matrix part = minor(m, r, c);
    result = pow(-1, (r + c + 2)) * determinant(&part); // + 2 is added as it is 0 based indexed
    return result;
}

Matrix cofactorMatrix(Matrix *m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Cofactors operations valid on Square Matrices only");
    }

    Matrix cofactor_matrix = {(double *)malloc(sizeof(double) * m->rows * m->cols), m->rows, m->cols};

    if (cofactor_matrix.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed for Cofactor Matrix");
    }

    for (unsigned int i = 0; i < m->rows; ++i)
    {
        for (unsigned int j = 0; j < m->cols; ++j)
        {
            double result = cofactor(m, i, j);
            *(cofactor_matrix.data + i * m->cols + j) = result;
        }
    }

    return cofactor_matrix;
}

Matrix adjoint(Matrix *m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Adjoint operations valid on Square Matrices only");
    }

    Matrix adjoint = {(double *)malloc(sizeof(double) * m->rows * m->cols), m->rows, m->cols};

    if (adjoint.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed for Adjoint Matrix");
    }

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            Matrix minorMatrix = minor(m, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            int det = determinant(&minorMatrix);
            *(adjoint.data + j * adjoint.cols + i) = sign * det;
            free(minorMatrix.data);
        }
    }

    return adjoint;
}

Matrix inverse(Matrix *m)
{
    if (!isSquareMatrix(m))
    {
        HANDLE_ERROR_MSG("Inverse operation valid on Square Matrices only");
    }

    double det = determinant(m);

    if (det == 0)
    {
        HANDLE_ERROR_MSG("Inverse does not exist for a singular matrix (determinant is zero)");
    }

    Matrix adj = adjoint(m);
    Matrix inverse_matrix = {(double *)malloc(sizeof(double) * m->rows * m->cols), m->rows, m->cols};

    if (inverse_matrix.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed for Inverse Matrix");
    }

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            *(inverse_matrix.data + i * inverse_matrix.cols + j) = *(adj.data + i * adj.cols + j) / (double)det;
        }
    }

    return inverse_matrix;
}

void test(void);

int main(void)
{
    srand(3);
    test();
    return 0;
}

void test(void)
{
    double m_[9];
    Matrix m = {m_, 3, 3};
    randomizeMatrix(&m, 100, 0);
    viewMatrix(&m, 1);

    getMatrixOrder(&m);
    printf("TOTAL: %d\n", getTotalElements(&m));
    printf("A(1, 1): %d\n", getMatrixElement(&m, 1, 1));
    printf("|A|: %f\n", determinant(&m));

    Matrix trp = transpose(&m);
    Matrix mMat = minorMatrix(&m);
    Matrix cMat = cofactorMatrix(&m);
    Matrix adj = adjoint(&m);
    Matrix inv = inverse(&m);

    viewMatrix(&trp, 1);
    viewMatrix(&mMat, 1);
    viewMatrix(&cMat, 1);
    viewMatrix(&adj, 1);
    viewMatrix(&inv, 1);
}