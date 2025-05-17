#include "Linlib.h"
using namespace Lin;

int main(void) {
    vec vector = {1, 2, 3, 4};
    printf("%d\n", vector[-1]);
    for (int i : vector) {
        printf("%d", i);
    }

    printf("\n");

    mat matrix = {{1, 2, 3}, {4, 5, 6}};
    printf("%d\n", matrix[-1][-1]);
    for (auto i : matrix) {
        printf("%d", i);
    }
}