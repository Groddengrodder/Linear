#include "Linlib.h"

int main(void) {
    auto vec1 = Lin::vec<int>(3).fill(2);
    auto vec2 = Lin::vec<int>(3).fill(3);
    auto matrix = Lin::mat<int>(3, 3);

    for (uint i = 0; i < matrix.getRows(); i++) {
        matrix[i][i] = 1;
    }
    vec1 += vec2;
    vec1 = matrix(vec1);
    vec1.print();
}
