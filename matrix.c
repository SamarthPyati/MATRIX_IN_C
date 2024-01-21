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

void getMatrixOrder(Matrix *matrix)
{
    printf("%zu x %zu\n", matrix->rows, matrix->cols);
}

size_t getTotalElements(Matrix *m)
{
    return m->rows * m->cols;
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
            MAT_ATP(matrix, i, j) = 0;
        }
    }
}

// row Major
void viewMatrix(Matrix *matrix, int addSpace, int viewAsINT)
{
    printf("[");
    for (int i = 0; i < matrix->rows; i++)
    {
        printf("[");
        for (int j = 0; j < matrix->cols; j++)
        {
            double el = MAT_ATP(matrix, i, j);

            // TODO: switch statement for different viewing types
            if (viewAsINT)
                addSpace ? printf(" %d ", (int)el) : printf("%d", (int)el);
            else
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

    Matrix result = mat_alloc(a->rows, a->cols);

    for (unsigned int i = 0; i < a->rows; i++)
    {
        for (unsigned int j = 0; j < a->cols; j++)
        {
            // *(result.data + i * result.cols + j) = *(a->data + i * a->cols + j) + *(b->data + i * b->cols + j);
            MAT_AT(result, i, j) = MAT_ATP(a, i, j) + MAT_ATP(b, i, j);
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

    Matrix result = mat_alloc(a->rows, a->cols);

    for (unsigned int i = 0; i < a->rows; i++)
    {
        for (unsigned int j = 0; j < a->cols; j++)
        {
            // *(result.data + i * result.cols + j) = *(a->data + i * a->cols + j) - *(b->data + i * b->cols + j);
            MAT_AT(result, i, j) = MAT_ATP(a, i, j) - MAT_ATP(b, i, j);
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

    Matrix result = mat_alloc(a->rows, b->cols);

    for (unsigned int i = 0; i < rows; i++)
    {
        for (unsigned int j = 0; j < cols; j++)
        {
            // sum = a[i][k] + b[k][j]
            int sum = 0;
            for (unsigned int k = 0; k < cols; k++)
            {
                // sum += (*(a->data + i * cols + k)) * (*(b->data + k * cols + j));
                sum += MAT_AT_EX_P(a, i, cols, k) * MAT_AT_EX_P(b, k, cols, j);
            }
            // *(result.data + i * cols + j) = sum;
            MAT_AT_EX(result, i, cols, j) = sum;
        }
    }
    return result;
}

Matrix transpose(Matrix *m)
{
    int rows = m->cols;
    int cols = m->rows;
    // A, A` | A[i][j] = A`[j][i]
    Matrix result = mat_alloc(rows, cols);
    // row & cols are swapped

    for (unsigned int i = 0; i < rows; ++i)
    {
        for (unsigned int j = 0; j < cols; ++j)
        {
            // *(result.data + i * cols + j) = *(m->data + j * rows + i);
            MAT_AT(result, i, j) = MAT_ATP(m, j, i);
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

    Matrix minor_ = mat_alloc(matrix->rows - 1, matrix->cols - 1);

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
            // *(minor_.data + minorRows * minor_.cols + minorCols) = *(matrix->data + i * matrix->cols + j);
            MAT_AT(minor_, minorRows, minorCols) = MAT_ATP(matrix, i, j);
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

    Matrix result = mat_alloc(m->rows, m->cols);

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            Matrix tempMinorMatrix = minor(m, i, j);
            double det = determinant(&tempMinorMatrix);
            // *(result.data + i * result.cols + j) = det;
            MAT_AT(result, i, j) = det;
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
            // *(cofactor_matrix.data + i * m->cols + j) = result;
            MAT_AT_EX(cofactor_matrix, i, m->cols, j) = result;
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

    Matrix adjoint = mat_alloc(m->rows, m->cols);

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            Matrix minorMatrix = minor(m, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            int det = determinant(&minorMatrix);
            // *(adjoint.data + j * adjoint.cols + i) = sign * det;
            MAT_AT(adjoint, j, i) = sign * det;
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
    Matrix inverse_matrix = mat_alloc(m->rows, m->cols);

    for (unsigned int i = 0; i < m->rows; i++)
    {
        for (unsigned int j = 0; j < m->cols; j++)
        {
            // *(inverse_matrix.data + i * inverse_matrix.cols + j) = *(adj.data + i * adj.cols + j) / (double)det;
            MAT_AT(inverse_matrix, i, j) = MAT_AT(adj, i, j) / det;
        }
    }

    return inverse_matrix;
}

Matrix identity(size_t rows, size_t cols)
{
    Matrix id = mat_alloc(rows, cols);
    clearMatrix(&id);
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            if (i == j)
            {
                MAT_AT(id, i, j) = 1;
            }
        }
    }
    return id;
}

void test(void);
void test_mat_at(void);

int main(void)
{
    srand(4); // Set seed value for pseudoRandom Generator
    test();
    return 0;
}

void test(void)
{
    /* FIRST WAY */
    Matrix m = mat_alloc(3, 3);
    randomizeMatrix(&m, 10, 1);
    viewMatrix(&m, 1, 1);

    /* GENERAL TESTING */
    // getMatrixOrder(&m);
    // printf("TOTAL: %d\n", getTotalElements(&m));
    // printf("A(1, 1): %d\n", (int)MAT_AT(m, 1, 1));
    // printf("|A|: %f\n\n", determinant(&m));

    Matrix trp = transpose(&m);
    Matrix mMat = minorMatrix(&m);
    Matrix cMat = cofactorMatrix(&m);
    Matrix adj = adjoint(&m);
    Matrix inv = inverse(&m);

    Matrix i = multiplyMatrix(&m, &inv);
    Matrix idt = identity(3, 3);
    viewMatrix(&trp, 1, 1);
    viewMatrix(&mMat, 1, 1);
    viewMatrix(&cMat, 1, 1);
    viewMatrix(&adj, 1, 1);
    viewMatrix(&inv, 1, 0);
    viewMatrix(&idt, 1, 1);
    viewMatrix(&i, 1, 1);
}

void test_mat_at(void)
{
    Matrix a = mat_alloc(3, 3);
    randomizeMatrix(&a, 1e1, 0);
    viewMatrix(&a, 1, 1);

    Matrix a_min = minorMatrix(&a);
    viewMatrix(&a_min, 1, 1);
    double a_ = determinant(&a_min);
    printf("DET(M(a)[0, 0]) : %d\n", (int)a_);
    // for (size_t i = 0; i < m.rows; ++i)
    // {
    //     for (size_t j = 0; j < m.cols; ++j)
    //     {
    //         printf("%d ", (int)MAT_ATP(&m, i, j));
    //     }
    //     printf("\n");
    // }
}