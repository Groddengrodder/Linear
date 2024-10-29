#include "Linlib.h"

const double epsilon = 1.e-10;
const uint dim = 10;

int main(void) {
    srand(time(NULL));

    Lin::mat<double> matrix(dim);
    while (matrix.det() < epsilon) {
        matrix.rand_fill(1.e8);
        matrix /= 1.e4;
    }

    Lin::mat<double> lu = Lin::LU(matrix);

    Lin::mat<double> u(dim);
    Lin::mat<double> l(dim);

    for (uint j = 0; j < dim; j++) {
        for (uint i = 0; i <= j; i++) {
            u[i][j] = lu[i][j];
        }

        for (uint i = j + 1; i < dim; i++) {
            l[i][j] = lu[i][j];
        }

        l[j][j] = 1;
    }

    Lin::mat<double> test = matrix - (l * u);
    double sum = 0;
    for (uint i = 0; i < dim * dim; i++) {
        sum += fabs(test[0][i]);
    }

    if (sum < epsilon) {
        printf("passed\n");
    } else {
        printf("failed\n");
    }

    Lin::mat<double> I(dim);
    while (I.det() < epsilon) {
        I.rand_fill(1.e8);
        I /= 1.e4;
    }

    auto b = Lin::vec<double>(dim).rand_fill(1.e8);
    b /= 1.e4;

    Lin::vec<double> r = I * Lin::solve_LU(I, b) - b;
    sum = 0;
    for (uint i = 0; i < r.size(); i++) {
        sum += fabs(r[i]);
    }

    if (sum < epsilon) {
        printf("passed\n");
    } else {
        printf("failed\n");
        I.print();
        printf("\n");
        Lin::solve_LU(I, b).print();
        printf("\n");
        r.print();
    }
}
