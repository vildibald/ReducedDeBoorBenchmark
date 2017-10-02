#pragma once

#include <vector>
#include "SurfaceDimension.h"
#include "MathFunction.h"


class Spline {
    SurfaceDimension xDimension_;
    SurfaceDimension yDimension_;
    std::vector<std::vector<double>> z;
    std::vector<std::vector<double>> dx;
    std::vector<std::vector<double>> dy;
    std::vector<std::vector<double>> dxy;

    Spline();

public:
    static Spline Null();

    bool IsNull();

    Spline(SurfaceDimension rowDimension, SurfaceDimension columnDimension);

    size_t RowsCount() const {
        return xDimension_.knotCount;
    }

    size_t ColumnsCount() const {
        return yDimension_.knotCount;
    }

    double X(const size_t j) const {
        auto h = abs(xDimension_.max - xDimension_.min)/(xDimension_.knotCount - 1);
        return j*h + xDimension_.min;
    }

    double Y(const size_t i) const {
        auto h = abs(yDimension_.max - yDimension_.min)/(yDimension_.knotCount - 1);
        return i*h + yDimension_.min;
    }

    const SurfaceDimension& XDimensionParameters() const {
        return xDimension_;
    }

    const SurfaceDimension& YDimensionParameters() const {
        return yDimension_;
    }

    double Z(const size_t i, const size_t j) const {
        return z[i][j];
    }

    double Dx(const size_t i, const size_t j) const {
        return dx[j][i];
    }

    double Dy(const size_t i, const size_t j) const {
        return dy[i][j];
    }

    double Dxy(const size_t i, const size_t j) const {
        return dxy[i][j];
    }

    void SetZ(const size_t i, const size_t j, const double value) {
        z[i][j] = value;
    }

    void SetDx(const size_t i, const size_t j, const double value) {
        dx[j][i] = value;
    }

    void SetDy(const size_t i, const size_t j, const double value) {
        dy[i][j] = value;
    }

    void SetDxy(const size_t i, const size_t j, const double value) {
        dxy[i][j] = value;
    }

    void Initialize(InterpolativeMathFunction mathFunction);

    void Print();
};

