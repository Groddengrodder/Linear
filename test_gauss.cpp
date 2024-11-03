#include "Linlib.h"
#include <time.h>

int main(int argc, char *argv[]) {
    uint dim = 1000;
    srand(time(NULL));
    const double epsilon = 1.e-7;

    if (argc == 2) {
        dim = atoi(argv[1]);
    }

    auto matrix = Lin::mat<double>(dim).rand_fill(100);
    auto vector = Lin::vec<double>(dim).rand_fill(100);

    auto s = clock();
    Lin::vec<double> solution = Lin::solve_gauss(matrix, vector);
    auto e = clock();

    Lin::vec<double> residuum = matrix * solution;
    residuum -= vector;

    double max = 0;
    for (uint i = 0; i < residuum.size(); i++) {
        if (fabs(residuum[i]) > max) {
            max = residuum[i];
        }
    }

    if (max < epsilon) {
        printf("passed\n");
    } else {
        printf("matrix:\n");
        matrix.print();
        printf("vector b:\n");
        vector.print();
        printf("solution:\n");
        solution.print();
        printf("residuum:\n");
        residuum.print();
        printf("max error = %.2le\n", max);
    }
    printf("time = %lf s\n", (double)(e - s) / CLOCKS_PER_SEC);
}
