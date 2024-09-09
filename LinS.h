#ifndef __LINS__
#define __LINS__
#include "Linlib.h"

namespace Lin {
class vec3 : public vec<double> {
    private:
    inline static const uint dim = 3;

    public:
    double &x;
    double &y;
    double &z;

    vec3() : vec<double>(dim), x(comp[0]), y(comp[1]), z(comp[2]) {}
    vec3(double input) : vec<double>(dim), x(comp[0]), y(comp[1]), z(comp[2]) {
        for (uint i = 0; i < dim; i++) {
            comp[i] = input;
        }
    }
    vec3(double d1, double d2, double d3) : vec<double>(dim), x(comp[0]), y(comp[1]), z(comp[2]) {
        comp[0] = d1;
        comp[1] = d2;
        comp[3] = d3;
    }
};

class vec2 : public vec<double> {
    private:
    inline static const uint dim = 2;

    public:
    double &x;
    double &y;

    vec2() : vec<double>(dim), x(comp[0]), y(comp[1]) {}
    vec2(double input) : vec<double>(dim), x(comp[0]), y(comp[1]) {
        for (uint i = 0; i < dim; i++) {
            comp[i] = input;
        }
    }
    vec2(double d1, double d2) : vec<double>(dim), x(comp[0]), y(comp[1]) {
        comp[0] = d1;
        comp[1] = d2;
    }
};

class mat3 : public mat<double> {
    private:
    inline static const uint dim = 3;

    public:
    mat3() : mat<double>(dim) {}
    mat3(double input) : mat<double>(dim) {
        for (uint i = 0; i < dim * dim; i++) {
            comp[0][i] = input;
        }
    }
};

class mat2 : public mat<double> {
    private:
    inline static const uint dim = 2;

    public:
    double &xx;
    double &xy;
    double &yx;
    double &yy;

    mat2() : mat<double>(dim), xx(comp[0][0]), xy(comp[0][1]), yx(comp[1][0]), yy(comp[1][1]) {}
    mat2(double input)
        : mat<double>(dim), xx(comp[0][0]), xy(comp[0][1]), yx(comp[1][0]), yy(comp[1][1]) {
        for (uint i = 0; i < dim * dim; i++) {
            comp[0][i] = input;
        }
    }
    mat2(double XX, double XY, double YX, double YY)
        : mat<double>(dim), xx(comp[0][0]), xy(comp[0][1]), yx(comp[1][0]), yy(comp[1][1]) {
        comp[0][0] = XX;
        comp[0][1] = XY;
        comp[1][0] = YX;
        comp[1][1] = YY;
    }
};
} // namespace Lin
#endif
