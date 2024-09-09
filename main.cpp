#include "LinS.h"

int main(void) {
    Lin::vec3 myVec;
    myVec.comp[0] = 0;
    myVec.comp[1] = 5;
    myVec.comp[2] = 4;

    myVec.x = 2;

    Lin::mat3 myMat(1);

    Lin::mat2 matrix = Lin::mat2(1, 0, 0, 1);

    print(matrix);
    printf("\n");
    print(myVec);
    printf("\n");
    print(myMat);
    printf("\n");
    print(myMat * myVec);
}
