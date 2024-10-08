#ifndef __LINLIB__
#define __LINLIB__

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

        size = input_size;
        memset(comp, 0, size * sizeof(type));
    }

    vec(const vec<type> &vector) {
        comp = (type *)malloc(vector.getSize() * sizeof(type));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        size = vector.getSize();
        for (uint i = 0; i < size; i++) {
            comp[i] = vector.comp[i];
        }
    }

    vec(vec<type> &&vector) {
        size = vector.size;
        comp = vector.comp;

        vector.size = 0;
        vector.comp = NULL;
    }

    vec cross(const vec &other) const;
    type len() const;
    void print() const;
    vec &fill(type input);
    type &operator[](uint index) const;

    vec &operator+=(const vec &other);
    vec &operator-=(const vec &other);
    vec &operator*=(type other);
    vec &operator/=(type other);

    vec &operator=(const vec &vector) {
        uint Size = size < vector.getSize() ? size : vector.getSize();

        for (uint i = 0; i < Size; i++) {
            comp[i] = vector.comp[i];
        }

        for (uint i = Size; i < size; i++) {
            comp[i] = 0;
        }

        return *this;
    }

    vec &operator=(vec &&vector) {
        if (this == &vector) {
            return *this;
        }

        size = vector.size;
        if (comp != NULL) {
            type *temp = comp;
            free(temp);
        }

        comp = vector.comp;

        vector.size = 0;
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

    uint getSize() const { return size; }

    private:
    uint size;
};

template <typename type> class mat {
    public:
    type **comp;

    mat(uint input_rows, uint input_columns) {
        columns = input_columns;
        rows = input_rows;

        comp = (type **)malloc(rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(rows * columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < rows; i++) {
            comp[i] = col + i * columns;
        }

        memset(col, 0, rows * columns * sizeof(type));
    }

    mat(uint input) {
        columns = input;
        rows = input;

        comp = (type **)malloc(rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(rows * columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < rows; i++) {
            comp[i] = col + i * columns;
        }

        memset(col, 0, rows * columns * sizeof(type));
    }

    mat(const mat<type> &matrix) {
        columns = matrix.getColumns();
        rows = matrix.getRows();

        comp = (type **)malloc(rows * sizeof(type *));
        if (comp == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        col = (type *)malloc(rows * columns * sizeof(type));
        if (col == NULL) {
            printf("Error: Couldnt allocate memory\n");
            exit(1);
        }

        for (uint i = 0; i < rows; i++) {
            comp[i] = col + i * columns;
        }

        for (uint i = 0; i < rows; i++) {
            for (uint j = 0; j < columns; j++) {
                comp[i][j] = matrix.comp[i][j];
            }
        }
    }

    mat(mat<type> &&matrix) {
        columns = matrix.columns;
        rows = matrix.rows;

        col = matrix.col;
        comp = matrix.comp;

        matrix.col = NULL;
        matrix.comp = NULL;
        matrix.columns = 0;
        matrix.rows = 0;
    }

    type tr() const;
    type det() const;
    void print() const;
    mat &fill(type input);
    type *operator[](uint index) const;
    vec<type> operator()(const vec<type> &vector) const;

    mat &operator+=(const mat &other);
    mat &operator-=(const mat &other);
    mat &operator*=(type other);
    mat &operator/=(type other);

    mat &operator=(const mat &matrix) {
        if (rows != matrix.rows || columns != matrix.columns) {
            printf("please\n");
            exit(1);
        }

        for (uint i = 0; i < rows; i++) {
            for (uint j = 0; j < columns; j++) {
                comp[i][j] = matrix.comp[i][j];
            }
        }

        return *this;
    }

    mat &operator=(mat &&matrix) {
        if (this == &matrix) {
            return *this;
        }

        rows = matrix.rows;
        columns = matrix.columns;

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

        matrix.rows = 0;
        matrix.columns = 0;
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

    uint getRows() const { return rows; }

    uint getColumns() const { return columns; }

    private:
    type *col;

    uint columns;
    uint rows;
};

template <typename type> void vec<type>::print() const {
    for (uint i = 0; i < size; i++) {
        printf("%.2le\n", (double)comp[i]);
    }
}

template <typename type> void mat<type>::print() const {
    for (uint i = 0; i < rows; i++) {
        for (uint j = 0; j < columns; j++) {
            printf("%.2le ", (double)comp[i][j]);
        }
        printf("\n");
    }
}

template <typename type> void is_same_size(const vec<type> &vec1, const vec<type> &vec2) {
    if (vec1.getSize() != vec2.getSize()) {
        fprintf(stderr, "Error: Vectors have to have the same size\n");
        exit(1);
    }
}

template <typename type> void is_same_size_mult(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.getColumns() != mat2.getRows()) {
        fprintf(stderr, "Error: Matrices have to have the same size\n");
        exit(1);
    }
}

template <typename type> void is_same_size_add(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.getRows() != mat2.getRows() || mat1.getColumns() != mat2.getColumns()) {
        fprintf(stderr, "Error: Matrices have to have the same size\n");
        exit(1);
    }
}

template <typename type> vec<type> &vec<type>::fill(type input) {
    for (uint i = 0; i < size; i++) {
        comp[i] = input;
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator+=(const vec<type> &other) {
    is_same_size(*this, other);

    for (uint i = 0; i < size; i++) {
        comp[i] += other[i];
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator-=(const vec<type> &other) {
    is_same_size(*this, other);

    for (uint i = 0; i < size; i++) {
        comp[i] -= other[i];
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator*=(type other) {
    for (uint i = 0; i < size; i++) {
        comp[i] *= other;
    }

    return *this;
}

template <typename type> vec<type> &vec<type>::operator/=(type other) {
    if (other == 0) {
        printf("dont divide by 0\n");
        exit(1);
    }

    for (uint i = 0; i < size; i++) {
        comp[i] /= other;
    }

    return *this;
}

template <typename type> bool operator==(const vec<type> &vec1, const vec<type> &vec2) {
    if (vec1.getSize() != vec2.getSize()) {
        return false;
    }

    for (uint i = 0; i < vec1.getSize(); i++) {
        if (vec1.comp[i] != vec2.comp[i]) {
            return false;
        }
    }

    return true;
}

template <typename type> vec<type> operator+(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    vec<type> solution(vec1.getSize());
    for (uint i = 0; i < vec1.getSize(); i++) {
        solution.comp[i] = vec1.comp[i] + vec2.comp[i];
    }

    return solution;
}

template <typename type> vec<type> operator-(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    vec<type> solution(vec1.getSize());
    for (uint i = 0; i < vec1.getSize(); i++) {
        solution.comp[i] = vec1.comp[i] - vec2.comp[i];
    }

    return solution;
}

template <typename type> vec<type> operator*(const type &scalar, const vec<type> &vector) {
    vec<type> solution(vector.getSize());
    for (uint i = 0; i < vector.getSize(); i++) {
        solution.comp[i] = scalar * vector.comp[i];
    }

    return solution;
}

template <typename type> type vec<type>::len() const {
    double solution = 0;

    for (uint i = 0; i < size; i++) {
        solution += pow(comp[i], 2);
    }
    solution = sqrt(solution);

    return (type)solution;
}

template <typename type> vec<type> vec<type>::cross(const vec<type> &other) const {
    if (size != 3 || other.size != 3) {
        printf("Error: Cross product only works in d=3\n");
    }

    vec<type> solution(3);

    solution.comp[0] = comp[1] * other.comp[2] - comp[2] * other.comp[1];
    solution.comp[1] = comp[2] * other.comp[0] - comp[0] * other.comp[2];
    solution.comp[2] = comp[0] * other.comp[1] - comp[1] * other.comp[0];

    return solution;
}

template <typename type> type &vec<type>::operator[](uint index) const {
    if (index >= size) {
        printf("Attempted to acess out of bounds element in vector\n");
        exit(1);
    }

    return comp[index];
}

template <typename type> type operator*(const vec<type> &vec1, const vec<type> &vec2) {
    is_same_size(vec1, vec2);

    type solution = 0;
    for (uint i = 0; i < vec1.getSize(); i++) {
        solution += vec1.comp[i] * vec2.comp[i];
    }

    return solution;
}

template <typename type> bool operator==(const mat<type> &mat1, const mat<type> &mat2) {
    if (mat1.getRows() != mat2.getRows() || mat1.getColumns() != mat2.getColumns()) {
        return false;
    }

    for (uint i = 0; i < mat1.getRows(); i++) {
        for (uint j = 0; j < mat1.getColumns(); j++) {
            if (mat1.comp[i][j] != mat2.comp[i][j]) {
                return false;
            }
        }
    }

    return true;
}

template <typename type> mat<type> &mat<type>::fill(type input) {
    for (uint i = 0; i < rows * columns; i++) {
        comp[0][i] = input;
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator+=(const mat<type> &other) {
    is_same_size_add(*this, other);

    for (uint i = 0; i < rows * columns; i++) {
        comp[0][i] += other[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator-=(const mat<type> &other) {
    is_same_size_add(*this, other);

    for (uint i = 0; i < rows * columns; i++) {
        comp[0][i] -= other[0][i];
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator*=(type other) {
    for (uint i = 0; i < rows * columns; i++) {
        comp[0][i] *= other;
    }

    return *this;
}

template <typename type> mat<type> &mat<type>::operator/=(type other) {
    for (uint i = 0; i < rows * columns; i++) {
        comp[0][i] /= other;
    }

    return *this;
}

template <typename type> mat<type> operator*(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_mult(mat1, mat2);

    mat<type> solution(mat1.getRows(), mat2.getColumns());
    for (uint i = 0; i < mat1.getRows(); i++) {
        for (uint j = 0; j < mat2.getColumns(); j++) {
            for (uint k = 0; k < mat1.getColumns(); k++) {
                solution.comp[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return solution;
}

template <typename type> mat<type> operator*(const type scalar, const mat<type> &matrix) {
    mat<type> solution(matrix.getRows(), matrix.getColumns());
    for (uint i = 0; i < matrix.getRows(); i++) {
        for (uint j = 0; j < matrix.getColumns(); j++) {
            solution.comp[i][j] = scalar * matrix.comp[i][j];
        }
    }

    return solution;
}

template <typename type> mat<type> operator+(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_add(mat1, mat2);

    mat<type> solution(mat1.getRows(), mat1.getColumns());
    for (uint i = 0; i < mat1.getRows(); i++) {
        for (uint j = 0; j < mat1.getColumns(); j++) {
            solution.comp[i][j] = mat1.comp[i][j] + mat2.comp[i][j];
        }
    }

    return solution;
}

template <typename type> mat<type> operator-(const mat<type> &mat1, const mat<type> &mat2) {
    is_same_size_add(mat1, mat2);

    mat<type> solution(mat1.getRows(), mat1.getColumns());
    for (uint i = 0; i < mat1.getRows(); i++) {
        for (uint j = 0; j < mat1.getColumns(); j++) {
            solution.comp[i][j] = mat1.comp[i][j] - mat2.comp[i][j];
        }
    }

    return solution;
}

template <typename type> vec<type> operator*(const mat<type> &matrix, const vec<type> &vector) {
    vec<type> solution(vector.getSize());
    if (vector.getSize() != matrix.getColumns()) {
        printf("welp\n");
        exit(1);
    }

    for (uint i = 0; i < matrix.getRows(); i++) {
        for (uint j = 0; j < matrix.getColumns(); j++) {
            solution.comp[i] += matrix.comp[i][j] * vector.comp[j];
        }
    }

    return solution;
}

template <typename type> vec<type> mat<type>::operator()(const vec<type> &vector) const {
    return ((*this) * vector);
}

template <typename type> type mat<type>::tr() const {
    if (rows != columns) {
        printf("Error: need square matrix\n");
        exit(1);
    }

    type solution = 0;

    for (uint i = 0; i < rows; i++) {
        solution += comp[i][i];
    }

    return solution;
}

template <typename type> type *mat<type>::operator[](uint index) const { return comp[index]; }

template <typename type> mat<type> mat_minor(const mat<type> &matrix, uint row, uint column) {
    mat<type> minor(matrix.getRows() - 1, matrix.getColumns() - 1);

    bool offi = 0;
    bool offj = 0;

    for (uint i = 0; i < minor.rows; i++) {
        if (i >= row) {
            offi = 1;
        }
        offj = 0;
        for (uint j = 0; j < minor.columns; j++) {
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

    if (rows != columns) {
        printf("Error: Determinant only defined on square matrices\n");
        exit(1);
    }

    if (rows == 1) {
        determinant = comp[0][0];
    } else {
        int one = 1;
        for (uint i = 0; i < rows; i++) {
            determinant += one * comp[i][0] * mat_minor(*this, i, 0).det();
            one *= -1;
        }
    }
    return determinant;
}

template <typename type> vec<type> vec_solve(const mat<type> &matrix, const vec<type> &vector) {
    mat<type> matrix_temp(matrix.getRows(), matrix.getColumns());
    vec<type> vector_temp(vector.getSize());
    vec<type> solution(vector.getSize());

    if (matrix.getColumns() != vector.rows) {
        printf("Error: Invalid system of equations\n");
        exit(1);
    }

    for (uint i = 0; i < matrix.getRows(); i++) {
        for (uint j = 0; j < matrix.getColumns(); j++) {
            matrix_temp.comp[i][j] = matrix.comp[i][j];
        }
    }
    for (uint i = 0; i < matrix.getRows(); i++) {
        vector_temp.comp[i] = vector.comp[i];
        solution.comp[i] = 0;
    }

    for (uint i = 0; i < matrix.getColumns(); i++) {
        int temp = -1;
        for (uint j = 0; j < matrix.getRows(); j++) {
            for (uint k = 0; k < matrix.getColumns(); k++) {
                if (matrix_temp.comp[j][k] != 0 && k != i) {
                    break;
                }
                if (matrix_temp.comp[j][k] != 0 && k == i) {
                    temp = j;
                    break;
                }
            }
            if (temp != -1) {
                break;
            }
        }
        if (temp != -1) {
            for (uint j = 0; j < matrix.getRows(); j++) {
                if (temp != j) {
                    if (matrix_temp.comp[j][i] != 0) {
                        for (uint k = 0; k < matrix.getColumns(); k++) {
                            matrix_temp.comp[j][k] =
                                matrix_temp.comp[j][k] + matrix_temp.comp[temp][k];
                        }
                        vector_temp.comp[j] = vector_temp.comp[j] + vector_temp.comp[temp];
                    }
                }
            }
        }
    }

    for (uint i = 0; i < matrix.getRows(); i++) {
        int temp = -2;
        for (uint j = 0; j < matrix.getColumns(); j++) {
            if (matrix_temp.comp[i][j] != 0) {
                temp = j;
                break;
            }
        }
        if (temp != -2) {
            solution.comp[temp] = vector_temp.comp[i];
        } else if (temp == -2 && vector_temp.comp[i] != 0) {
            printf("Warning: System of equations was unsolvable\n");
            break;
        }
    }

    return solution;
}
} // namespace Lin
#endif
