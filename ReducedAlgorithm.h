#pragma once

#include "MathFunction.h"
#include "Spline.h"
#include "SurfaceDimension.h"
#include "Tridiagonal.h"
#include <functional>
#include <vector>
#include <omp.h>

class ReducedAlgorithm {
    SurfaceDimension xDimension_, yDimension_;
    InterpolativeMathFunction function_;
    Tridiagonals xTridiagonals_;
    Tridiagonals yTridiagonals_;
    Timer timer_;
    bool isParallel_;

public:
    ReducedAlgorithm(const SurfaceDimension xDimension,
                     const SurfaceDimension yDimension,
                     const InterpolativeMathFunction f);

    Spline Calculate();

    double ExecutionTime();

    double AllTime();

    void InParallel(bool value);

    bool IsParallel();

private:
    void Initialize(Spline& values);

    void FillDx(Spline& values);

    void FillDy(Spline& values);

    void FillDxy(Spline& values);

    void FillDyx(Spline& values);

    void Parallelize(bool inParallel);

    void InitializeTridiagonals(Spline& spline);
};


