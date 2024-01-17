/**
    Matrix Multiplication
    matrices.c
    Matrix data structure in C.

    @author Samarth Pyati
    @version 1.0 16/1/24
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
    long int *data;
    int rows;
    int cols;
} Matrix;

void getMatrixOrder(Matrix *matrix)
{
    printf("%d x %d\n", matrix->rows, matrix->cols);
}

int validateMatrix(Matrix matrix)
{
    return 0;
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
            long int el = *(matrix->data + i * matrix->cols + j);
            addSpace ? printf(" %ld ", el) : printf("%ld", el);

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

    Matrix result = {(long int *)malloc(sizeof(long int) * a->rows * a->cols), a->rows, a->cols};

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

    Matrix result = {(long int *)malloc(sizeof(long int) * a->rows * a->cols), a->rows, a->cols};

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
    Matrix result = {malloc(sizeof(long int) * a->rows * b->cols), a->rows, b->cols};

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
    Matrix result = {malloc(sizeof(long int) * rows * cols), rows, cols};
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
Matrix minorMatrix(Matrix *matrix, int r, int c)
{
    if (r < 0 || r >= matrix->rows && c < 0 || c >= matrix->cols)
    {
        HANDLE_ERROR_MSG("Given row or column are out of bounds of the dimensions of Matrix.");
    }

    Matrix minor = {malloc(sizeof(long int) * (matrix->rows - 1) * (matrix->cols - 1)), matrix->rows - 1, matrix->cols - 1};

    if (minor.data == NULL)
    {
        HANDLE_ERROR_MSG("Memory Allocation Failed For Matrix");
    }

    int minorRows = 0, minorCols = 0;

    for (int i = 0; i < matrix->rows; i++)
    {
        if (i == r)
        {
            continue; // skipping the row of minor element
        }

        minorCols = 0; // Reset minor column for each row in the original matrix

        for (int j = 0; j < matrix->cols; j++)
        {
            if (j == c)
            {
                continue;
            }

            // copying the desired data
            *(minor.data + minorRows * minor.cols + minorCols) = *(matrix->data + i * matrix->cols + j);
            minorCols++;
        }
        minorRows++;
    }
    return minor;
}

int main(void)
{
    srand(time(NULL));

    long int m_[9];

    Matrix m = {m_, 3, 3};
    randomizeMatrix(&m, 100, 10);
    Matrix a11 = minorMatrix(&m, 0, 0);
    viewMatrix(&m, 1);
    viewMatrix(&a11, 1);
    return 0;
}