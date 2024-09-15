#include "Linlib.h"

int main(void) {
    Lin::vec<int> myVec(5);
    int a = myVec[3];
    myVec.cross(Lin::vec<int>(5));
    Lin::mat<int> matrix(4, 4);
    matrix[1][1];
}
