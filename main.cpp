#include "Linlib.h"

int main(void) {
    Lin::mat<double> matrix(3);
    matrix[0][0] = 1;
    for (uint i = 1; i < 9; i++) {
        matrix[0][i] = i;
    }

    Lin::mat<double> lu = Lin::LU(matrix);

    Lin::mat<double> u(3);
    Lin::mat<double> l(3);

    for (uint j = 0; j < 3; j++) {
        for (uint i = 0; i <= j; i++) {
            u[i][j] = lu[i][j];
        }

        for (uint i = j + 1; i < 3; i++) {
            l[i][j] = lu[i][j];
        }

        l[j][j] = 1;
    }

    Lin::mat<double> test = matrix - (l * u);
    double sum = 0;
    for (uint i = 0; i < 9; i++) {
        sum += test[0][i];
    }

    if (!sum) {
        printf("passed\n");
    } else {
        printf("failed\n");
    }
}
