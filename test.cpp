#define __TRACK_ALLOC__
#include "Linlib.h"
using namespace Lin;

int main(void) {
    vec vector = {1, 2, 3, 4};
    printf("%d\n", vector[-1]);
    for (int i : vector) {
        printf("%d", i);
    }

    printf("\n");

    mat mat1 = {{1, 2, 3}, {4, 5, 6}};
    mat mat2 = {{1, 2, 3}, {4, 5, 6}};
    printf("%d\n", mat1[-1][-1]);
    for (auto i : mat1) {
        printf("%d", i);
    }
    printf("\n");

    mat mat3 = mat1 + mat2;
    mat3.print();

    printf("allocs = %d\n", alloc_tracker);
}