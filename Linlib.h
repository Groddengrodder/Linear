#ifndef __LINLIB__
#define __LINLIB__

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace Lin {
template <typename type> class mat;

template <typename type> class vec {
    public:
    type *comp;

    vec(uint input_size) {
        comp = (type *)malloc(input_size * sizeof(type));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        m_size = input_size;
        memset(comp, 0, m_size * sizeof(type));
    }

    vec(const vec<type> &vector) {
        comp = (type *)malloc(vector.size() * sizeof(type));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        m_size = vector.size();
        for (uint i = 0; i < m_size; i++) {
            comp[i] = vector.comp[i];
        }
    }

    vec(vec<type> &&vector) {
        m_size = vector.m_size;
        comp = vector.comp;

        vector.m_size = 0;
        vector.comp = NULL;
    }

    vec cross(const vec &other) const;
    type len() const;
    void print() const;
    vec &fill(type input);
    vec &rand_fill(uint input);
    type &operator[](uint index) const;

    vec &operator+=(const vec &other);
    vec &operator-=(const vec &other);
    vec &operator*=(type other);
    vec &operator/=(type other);

    vec &add(const vec &vec1, const vec &vec2);
    vec &mult(type scalar, const vec &vector);
    vec &mult(const mat<type> &matrix, const vec &vector);

    vec &operator=(const vec &vector) {
        uint Size = m_size < vector.size() ? m_size : vector.size();

        for (uint i = 0; i < Size; i++) {
            comp[i] = vector.comp[i];
        }

        for (uint i = Size; i < m_size; i++) {
            comp[i] = 0;
        }

        return *this;
    }

    vec &operator=(vec &&vector) {
        if (this == &vector) {
            return *this;
        }

        m_size = vector.m_size;
        if (comp != NULL) {
            type *temp = comp;
            free(temp);
        }

        comp = vector.comp;

        vector.m_size = 0;
        vector.comp = NULL;

        return *this;
    }

    ~vec() {
        type *temp = comp;
        comp = NULL;
        if (temp != NULL) {
            free(temp);
        }
    }

    uint size() const { return m_size; }

    private:
    uint m_size;
};

template <typename type> class mat {
    public:
    type **comp;

    mat(uint input_rows, uint input_columns) {
        m_columns = input_columns;
        m_rows = input_rows;

        comp = (type **)malloc(m_rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(m_rows * m_columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < m_rows; i++) {
            comp[i] = col + i * m_columns;
        }

        memset(col, 0, m_rows * m_columns * sizeof(type));
    }

    mat(uint input) {
        m_columns = input;
        m_rows = input;

        comp = (type **)malloc(m_rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(m_rows * m_columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < m_rows; i++) {
            comp[i] = col + i * m_columns;
        }

        memset(col, 0, m_rows * m_columns * sizeof(type));
    }

    mat(const mat<type> &matrix) {
        m_columns = matrix.columns();
        m_rows = matrix.rows();

        comp = (type **)malloc(m_rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(m_rows * m_columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < m_rows; i++) {
            comp[i] = col + i * m_columns;
        }

        for (uint i = 0; i < m_rows; i++) {
            for (uint j = 0; j < m_columns; j++) {
                comp[i][j] = matrix.comp[i][j];
            }
        }
    }

    mat(mat<type> &&matrix) {
        m_columns = matrix.m_columns;
        m_rows = matrix.m_rows;

        col = matrix.col;
        comp = matrix.comp;

        matrix.col = NULL;
        matrix.comp = NULL;
        matrix.m_columns = 0;
        matrix.m_rows = 0;
    }

    type tr() const;
    type det() const;
    type det_LU() const;
    void print() const;
    mat &fill(type input);
    mat &rand_fill(uint input);
    static mat identity(uint size);
    mat &transpose();
    type *&operator[](uint index) const;
    vec<type> operator()(const vec<type> &vector) const;

    mat &operator+=(const mat &other);
    mat &operator-=(const mat &other);
    mat &operator*=(type other);
    mat &operator/=(type other);

    mat &add(const mat &mat1, const mat &mat2);
    mat &mult(type scalar, const mat &matrix);
    mat &mult(const mat &mat1, const mat &mat2);

    mat &operator=(const mat &matrix) {
        if (m_rows != matrix.m_rows || m_columns != matrix.m_columns) {
            printf("please\n");
            exit(1);
        }

        for (uint i = 0; i < m_rows; i++) {
            for (uint j = 0; j < m_columns; j++) {
                comp[i][j] = matrix.comp[i][j];
            }
        }

        return *this;
    }

    mat &operator=(mat &&matrix) {
        if (this == &matrix) {
            return *this;
        }

        m_rows = matrix.m_rows;
        m_columns = matrix.m_columns;

        col = matrix.col;

        if (comp[0] != NULL) {
            type *temp = comp[0];
            free(temp);
        }

        if (comp != NULL) {
            type **temp = comp;
            free(temp);
        }

        comp = matrix.comp;

        matrix.m_rows = 0;
        matrix.m_columns = 0;
        matrix.col = NULL;
        matrix.comp = NULL;

        return *this;
    }

    ~mat() {
        type **temp1 = comp;
        type *temp2 = col;
        col = NULL;
        comp = NULL;
        if (temp1 != NULL) {
            free(temp1);
        }
        if (temp2 != NULL) {
            free(temp2);
        }
    }

    uint rows() const { return m_rows; }

    uint columns() const { return m_columns; }

    private:
    type *col;

    uint m_columns;
    uint m_rows;
};

template <> void vec<int>::print() const {
    for (uint i = 0; i < m_size; i++) {
        printf("%+d\n", comp[i]);
    }
}

template <> void mat<int>::print() const {
    for (uint i = 0; i < m_rows; i++) {
        for (uint j = 0; j < m_columns; j++) {
            printf("%+d ", comp[i][j]);
        }
        printf("\n");
    }
}

template <typename type> void vec<type>::print() const {
    for (uint i = 0; i < m_size; i++) {
        printf("%+.2le\n", (double)comp[i]);
    }
}

template <typename type> void mat<type>::print() const {
    for (uint i = 0; i < m_rows; i++) {
        for (uint j = 0; j < m_columns; j++) {
            printf("%+.2le ", (double)comp[i][j]);
        }
        printf("\n");
    }
}

template <typename type> void is_same_size(const vec<type> &vec1, const vec<type> &vec2) {
    if (vec1.size() != vec2.size()) {
        fprintf(stderr, "Error: Vectors have to have the same size\n");
        exit(1);
    }
}

template <typename type> void is_same_size_mult(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.columns() != mat2.rows()) {
        fprintf(stderr, "Error: Matrices have to have the same size\n");
        exit(1);
    }
}

template <typename type> void is_same_size_add(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.rows() != mat2.rows() || mat1.columns() != mat2.columns()) {
        fprintf(stderr, "Error: Matrices have to have the same size\n");
        exit(1);
    }
}

template <typename type> vec<type> &vec<type>::add(const vec<type> &vec1, const vec<type> &vec2) {
    assert(vec1.size() == vec2.size() && m_size == vec1.size());

    for (uint i = 0; i < m_size; i++) {
        comp[i] = vec1[i] + vec2[i];
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::mult(type scalar, const vec<type> &vector) {
    assert(m_size == vector.size());

    for (uint i = 0; i < m_size; i++) {
        comp[i] = scalar * vector[i];
    }

    return *this;
}

template <typename type>
vec<type> &vec<type>::mult(const mat<type> &matrix, const vec<type> &vector) {
    assert(m_size == matrix.rows());
    assert(matrix.columns() == vector.size());

    for (uint i = 0; i < m_size; i++) {
        type sum = 0;
        for (uint j = 0; j < vector.size(); j++) {
            sum += matrix[i][j] * vector[j];
        }

        comp[i] = sum;
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::add(const mat<type> &mat1, const mat<type> &mat2) {
    assert(mat1.rows() == mat2.rows() && m_rows == mat1.rows());
    assert(mat1.columns() == mat2.columns && m_columns == mat1.columns());

    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] = mat1[0][i] + mat2[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::mult(const type scalar, const mat<type> &matrix) {
    assert(m_rows == matrix.rows() && m_columns == matrix.columns());

    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] = scalar * matrix[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::mult(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_mult(mat1, mat2);
    assert(m_rows == mat1.rows() && m_columns == mat2.columns());

    for (uint i = 0; i < mat1.rows(); i++) {
        for (uint j = 0; j < mat2.columns(); j++) {
            type sum = 0;
            for (uint k = 0; k < mat1.columns(); k++) {
                sum += mat1[i][k] * mat2[k][j];
            }
            comp[i][j] = sum;
        }
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::transpose() {
    if (m_rows == m_columns) {
        for (uint i = 0; i < m_rows; i++) {
            for (uint j = 0; j < i; j++) {
                type temp = comp[i][j];
                comp[i][j] = comp[j][i];
                comp[j][i] = temp;
            }
        }
        return *this;
    }

    mat<type> matrix(m_columns, m_rows);
    for (uint i = 0; i < matrix.rows(); i++) {
        for (uint j = 0; j < matrix.columns(); j++) {
            matrix[i][j] = comp[j][i];
        }
    }

    *this = std::move(matrix);

    return *this;
}

template <typename type> vec<type> &vec<type>::fill(type input) {
    for (uint i = 0; i < m_size; i++) {
        comp[i] = input;
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::rand_fill(uint input) {
    for (uint i = 0; i < m_size; i++) {
        comp[i] = (type)(rand() % input);
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator+=(const vec<type> &other) {
    is_same_size(*this, other);

    for (uint i = 0; i < m_size; i++) {
        comp[i] += other[i];
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator-=(const vec<type> &other) {
    is_same_size(*this, other);

    for (uint i = 0; i < m_size; i++) {
        comp[i] -= other[i];
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator*=(type other) {
    for (uint i = 0; i < m_size; i++) {
        comp[i] *= other;
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator/=(type other) {
    if (other == 0) {
        printf("dont divide by 0\n");
        exit(1);
    }

    for (uint i = 0; i < m_size; i++) {
        comp[i] /= other;
    }

    return *this;
}

template <typename type> bool operator==(const vec<type> &vec1, const vec<type> &vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }

    for (uint i = 0; i < vec1.size(); i++) {
        if (vec1.comp[i] != vec2.comp[i]) {
            return false;
        }
    }

    return true;
}

template <typename type> vec<type> operator+(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    vec<type> solution(vec1.size());
    for (uint i = 0; i < vec1.size(); i++) {
        solution.comp[i] = vec1.comp[i] + vec2.comp[i];
    }

    return solution;
}

template <typename type> vec<type> operator-(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    vec<type> solution(vec1.size());
    for (uint i = 0; i < vec1.size(); i++) {
        solution.comp[i] = vec1.comp[i] - vec2.comp[i];
    }

    return solution;
}

template <typename type> vec<type> operator*(const type &scalar, const vec<type> &vector) {
    vec<type> solution(vector.size());
    for (uint i = 0; i < vector.size(); i++) {
        solution.comp[i] = scalar * vector.comp[i];
    }

    return solution;
}

template <typename type> type vec<type>::len() const {
    double solution = 0;

    for (uint i = 0; i < m_size; i++) {
        solution += pow(comp[i], 2);
    }
    solution = sqrt(solution);

    return (type)solution;
}

template <typename type> vec<type> vec<type>::cross(const vec<type> &other) const {
    if (m_size != 3 || other.m_size != 3) {
        printf("Error: Cross product only works in d=3\n");
    }

    vec<type> solution(3);

    solution.comp[0] = comp[1] * other.comp[2] - comp[2] * other.comp[1];
    solution.comp[1] = comp[2] * other.comp[0] - comp[0] * other.comp[2];
    solution.comp[2] = comp[0] * other.comp[1] - comp[1] * other.comp[0];

    return solution;
}

template <typename type> type &vec<type>::operator[](uint index) const {
    if (index >= m_size) {
        printf("Attempted to acess out of bounds element in vector\n");
        exit(1);
    }

    return comp[index];
}

template <typename type> type operator*(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    type solution = 0;
    for (uint i = 0; i < vec1.size(); i++) {
        solution += vec1.comp[i] * vec2.comp[i];
    }

    return solution;
}

template <typename type> bool operator==(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.rows() != mat2.rows() || mat1.columns() != mat2.columns()) {
        return false;
    }

    for (uint i = 0; i < mat1.rows(); i++) {
        for (uint j = 0; j < mat1.columns(); j++) {
            if (mat1.comp[i][j] != mat2.comp[i][j]) {
                return false;
            }
        }
    }

    return true;
}

template <typename type> mat<type> &mat<type>::fill(type input) {
    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] = input;
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::rand_fill(uint input) {
    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] = (type)(rand() % input);
    }

    return *this;
}

template <typename type> mat<type> mat<type>::identity(uint size) {
    mat<type> I(size);
    for (uint i = 0; i < I.rows(); i++) {
        I[i][i] = 1;
    }

    return I;
}

template <typename type> mat<type> &mat<type>::operator+=(const mat<type> &other) {
    is_same_size_add(*this, other);

    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] += other[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator-=(const mat<type> &other) {
    is_same_size_add(*this, other);

    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] -= other[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator*=(type other) {
    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] *= other;
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator/=(type other) {
    for (uint i = 0; i < m_rows * m_columns; i++) {
        comp[0][i] /= other;
    }

    return *this;
}

template <typename type> mat<type> operator*(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_mult(mat1, mat2);

    mat<type> solution(mat1.rows(), mat2.columns());
    for (uint i = 0; i < mat1.rows(); i++) {
        for (uint j = 0; j < mat2.columns(); j++) {
            for (uint k = 0; k < mat1.columns(); k++) {
                solution.comp[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return solution;
}

template <typename type> mat<type> operator*(const type scalar, const mat<type> &matrix) {
    mat<type> solution(matrix.rows(), matrix.columns());
    for (uint i = 0; i < matrix.rows(); i++) {
        for (uint j = 0; j < matrix.columns(); j++) {
            solution.comp[i][j] = scalar * matrix.comp[i][j];
        }
    }

    return solution;
}

template <typename type> mat<type> operator+(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_add(mat1, mat2);

    mat<type> solution(mat1.rows(), mat1.columns());
    for (uint i = 0; i < mat1.rows(); i++) {
        for (uint j = 0; j < mat1.columns(); j++) {
            solution.comp[i][j] = mat1.comp[i][j] + mat2.comp[i][j];
        }
    }

    return solution;
}

template <typename type> mat<type> operator-(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_add(mat1, mat2);

    mat<type> solution(mat1.rows(), mat1.columns());
    for (uint i = 0; i < mat1.rows(); i++) {
        for (uint j = 0; j < mat1.columns(); j++) {
            solution.comp[i][j] = mat1.comp[i][j] - mat2.comp[i][j];
        }
    }

    return solution;
}

template <typename type> vec<type> operator*(const mat<type> &matrix, const vec<type> &vector) {
    vec<type> solution(vector.size());
    if (vector.size() != matrix.columns()) {
        printf("welp\n");
        exit(1);
    }

    for (uint i = 0; i < matrix.rows(); i++) {
        for (uint j = 0; j < matrix.columns(); j++) {
            solution.comp[i] += matrix.comp[i][j] * vector.comp[j];
        }
    }

    return solution;
}

template <typename type> vec<type> mat<type>::operator()(const vec<type> &vector) const {
    return ((*this) * vector);
}

template <typename type> type mat<type>::tr() const {
    if (m_rows != m_columns) {
        printf("Error: need square matrix\n");
        exit(1);
    }

    type solution = 0;

    for (uint i = 0; i < m_rows; i++) {
        solution += comp[i][i];
    }

    return solution;
}

template <typename type> type *&mat<type>::operator[](uint index) const { return comp[index]; }

template <typename type> mat<type> LU(const mat<type> &matrix) {
    assert(matrix.rows() == matrix.columns());
    mat<type> lu(matrix.rows(), matrix.columns());
    assert(matrix[0][0] != 0);
    for (uint j = 0; j < matrix.rows(); j++) {
        type sum = 0;

        for (uint i = 0; i < j; i++) {
            sum = 0;
            for (uint k = 0; k < i; k++) {
                sum += lu[i][k] * lu[k][j];
            }

            lu[i][j] = matrix[i][j] - sum;
        }

        sum = 0;
        for (uint k = 0; k < j; k++) {
            sum += lu[j][k] * lu[k][j];
        }

        lu[j][j] = matrix[j][j] - sum;

        for (uint i = j + 1; i < matrix.rows(); i++) {
            sum = 0;
            for (uint k = 0; k < i; k++) {
                sum += lu[i][k] * lu[k][j];
            }

            assert(lu[j][j] != 0);
            lu[i][j] = 1. / lu[j][j] * (matrix[i][j] - sum);
        }
    }

    return lu;
}

template <typename type> mat<type> mat_minor(const mat<type> &matrix, uint row, uint column) {
    mat<type> minor(matrix.rows() - 1, matrix.columns() - 1);

    bool offi = 0;
    bool offj = 0;

    for (uint i = 0; i < minor.rows(); i++) {
        if (i >= row) {
            offi = 1;
        }
        offj = 0;
        for (uint j = 0; j < minor.columns(); j++) {
            if (j >= column) {
                offj = 1;
            }
            minor.comp[i][j] = matrix.comp[i + offi][j + offj];
        }
    }

    return minor;
}

template <typename type> type mat<type>::det() const {
    type determinant = 0;

    if (m_rows != m_columns) {
        printf("Error: Determinant only defined on square matrices\n");
        exit(1);
    }

    if (m_rows == 1) {
        determinant = comp[0][0];
    } else {
        double one = 1;
        for (uint i = 0; i < m_rows; i++) {
            determinant += one * comp[i][0] * mat_minor(*this, i, 0).det();
            one *= -1;
        }
    }
    return determinant;
}

template <typename type> vec<type> solve_cramer(const mat<type> &matrix, const vec<type> &vector) {
    assert(matrix.rows() == matrix.columns());
    assert(matrix.columns() == vector.size());

    vec<type> solution(vector.size());

    type det = matrix.det();
    assert(det != 0);
    for (uint j = 0; j < solution.size(); j++) {
        type sum = 0;
        for (uint i = 0; i < solution.size(); i++) {
            sum += vector[i] * pow(-1, i + j) * mat_minor(matrix, i, j).det();
        }

        solution[j] = sum / det;
    }

    return solution;
}

template <typename type> vec<type> solve_LU(const mat<type> &matrix, const vec<type> &vector) {
    assert(matrix.rows() == matrix.columns());
    assert(matrix.columns() == vector.size());

    vec<type> solution(vector.size());
    vec<type> y(vector.size());
    mat<type> lu = LU(matrix);

    for (uint i = 0; i < y.size(); i++) {
        type sum = 0;
        for (uint j = 0; j < i; j++) {
            sum += lu[i][j] * y[j];
        }

        y[i] = vector[i] - sum;
    }

    for (int i = vector.size() - 1; i >= 0; i--) {
        type sum = 0;
        for (uint j = i + 1; j < vector.size(); j++) {
            sum += lu[i][j] * solution[j];
        }
        assert(lu[i][i] != 0);
        solution[i] = (y[i] - sum) / lu[i][i];
    }

    return solution;
}

template <typename type> vec<type> solve_gauss(const mat<type> &matrix, const vec<type> &vector) {
    mat<type> matrix_temp = matrix;
    vec<type> vector_temp = vector;
    vec<type> solution(vector.size());

    if (matrix.columns() != vector.size()) {
        printf("Error: Invalid system of equations\n");
        exit(1);
    }

    const double epsilon = 1.e-8;

    for (uint i = 0; i < matrix.columns(); i++) {
        int temp = -1;
        for (uint j = 0; j < matrix.rows(); j++) {
            for (uint k = 0; k < matrix.columns(); k++) {
                if (fabs(matrix_temp[j][k]) >= epsilon && k != i) {
                    break;
                }
                if (fabs(matrix_temp[j][k]) >= epsilon && k == i) {
                    temp = j;
                    break;
                }
            }
            if (temp != -1) {
                break;
            }
        }

        if (temp == -1) {
            continue;
        }

        for (uint j = 0; j < matrix.rows(); j++) {
            if (temp == j) {
                continue;
            }
            if (fabs(matrix_temp[j][i]) < epsilon) {
                continue;
            }

            type factor = matrix_temp[j][i];

            for (uint k = 0; k < matrix.columns(); k++) {
                matrix_temp[j][k] =
                    matrix_temp[j][k] - matrix_temp[temp][k] / matrix_temp[temp][i] * factor;
            }
            vector_temp[j] = vector_temp[j] - vector_temp[temp] / matrix_temp[temp][i] * factor;
        }
    }

    for (uint i = 0; i < matrix.rows(); i++) {
        int temp = -2;
        for (uint j = 0; j < matrix.columns(); j++) {
            if (fabs(matrix_temp[i][j]) >= epsilon) {
                temp = j;
                break;
            }
        }
        if (temp != -2) {
            solution[temp] = vector_temp[i] / matrix_temp[i][temp];
        } else if (temp == -2 && fabs(vector_temp[i]) >= epsilon) {
            printf("Warning: System of equations was unsolvable\n");
            break;
        }
    }

    return solution;
}

template <typename type> type mat<type>::det_LU() const {
    assert(m_rows == m_columns);

    mat<type> lu = LU(*this);

    type det = 1;
    for (uint i = 0; i < m_rows; i++) {
        det *= lu[i][i];
    }

    return det;
}
} // namespace Lin
#endif
