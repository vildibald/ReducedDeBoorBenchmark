#include "ReducedAlgorithm.h"
#include "StopWatch.h"
#include "utils.h"

ReducedAlgorithm::ReducedAlgorithm(const SurfaceDimension xDimension,
                                   const SurfaceDimension yDimension,
                                   const InterpolativeMathFunction f)
        : xDimension_(xDimension), yDimension_(yDimension), xTridiagonals_(), yTridiagonals_(),
          function_(f), isParallel_(false) {
}

Spline ReducedAlgorithm::Calculate() {
    Spline spline{xDimension_, yDimension_};
    Initialize(spline);

    timer_.Reset();
    timer_.Start();
    FillDx(spline);
    FillDy(spline);
    FillDxy(spline);
    FillDyx(spline);
    timer_.Stop();

    return spline;
}

void ReducedAlgorithm::FillDx(Spline& spline) {
    double threeDivH = 3/xDimension_.H();
    double twelveDivH = 12/xDimension_.H();
    auto unknownsCount = xDimension_.KnotCount();
    auto even = unknownsCount % 2 == 0;
    auto tau = even ? 0 : 2;
    auto eta = even ? -4 : 1;
    auto upsilon = even ? unknownsCount - 2 : unknownsCount - 3;
    auto equationsCount = even?  unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
    utils::For(0, static_cast<int>(spline.ColumnsCount()), 1, isParallel_, [&](int j) {
//    for (size_t j = 0; j < spline.ColumnsCount(); ++j) {
        auto& tridiagonal = xTridiagonals_.Get();
        auto& rightSideX = tridiagonal.RightSideBuffer();
        for (size_t i = 0; i < equationsCount; ++i) {
            auto i21 = 2*(i + 1);
            rightSideX[i] = threeDivH*(spline.Z(i21 + 2, j) - spline.Z(i21 - 2, j))
                            - twelveDivH*(spline.Z(i21 + 1, j) - spline.Z(i21 - 1, j));
        }
        rightSideX[0] -= spline.Dx(0, j);
        rightSideX[rightSideX.size() - 1] =
                threeDivH*(spline.Z(upsilon + tau, j) - spline.Z(upsilon - 2, j))
                - twelveDivH*(spline.Z(upsilon + 1, j) - spline.Z(upsilon - 1, j))
                -eta*spline.Dx(spline.RowsCount() - 1, j);

        tridiagonal.Solve();

        for (int i = 0; i < rightSideX.size(); ++i) {
            auto i21 = 2*(i + 1);
            spline.SetDx(i21, j, rightSideX[i]);
        }

        for (size_t i = 1; i < spline.RowsCount(); i += 2) {
            spline.SetDx(i, j,
                         0.25 * (threeDivH * (spline.Z(i + 1, j) -
                                                spline.Z(i - 1, j)) -
                                 spline.Dx(i + 1, j) -
                                 spline.Dx(i - 1, j))
            );
        }
    });
}

void ReducedAlgorithm::FillDy(Spline& spline) {
    double threeDivH = 3/yDimension_.H();
    double twelveDivH = 12/yDimension_.H();
    auto unknownsCount = yDimension_.KnotCount();
    auto even = unknownsCount % 2 == 0;
    auto tau = even ? 0 : 2;
    auto eta = even ? -4 : 1;
    auto upsilon = even ? unknownsCount : unknownsCount - 1;
    auto equationsCount = even ? unknownsCount / 2 - 1
                               : unknownsCount / 2;
    utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i) {
//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        auto& tridiagonal = yTridiagonals_.Get();
        auto& rightSideY = tridiagonal.RightSideBuffer();
        for (size_t j = 0; j < equationsCount-1; ++j) {
            auto j21 = 2*(j + 1);
            rightSideY[j] = threeDivH*(spline.Z(i, j21 + 2) - spline.Z(i, j21 - 2))
                            - twelveDivH*(spline.Z(i, j21 + 1) - spline.Z(i, j21 - 1));
        }
        rightSideY[0] -= spline.Dy(i, 0);
        rightSideY[rightSideY.size() - 1] =
                threeDivH*(spline.Z(i, upsilon + tau) - spline.Z(i, upsilon - 2))
                - twelveDivH*(spline.Z(i, upsilon + 1) - spline.Z(i, upsilon - 1))
                -eta*spline.Dy(i, spline.ColumnsCount() - 1);

        tridiagonal.Solve();

        for (int j = 0; j < rightSideY.size(); ++j) {
            auto j21 = 2*(j + 1);
            spline.SetDy(i, j21, rightSideY[j]);
        }

        for (size_t j = 1; j < spline.RowsCount(); j += 2) {
            spline.SetDy(i, j,
                         0.25 * (threeDivH * (spline.Z(i, j + 1) -
                                              spline.Z(i, j - 1)) -
                                 spline.Dy(i, j + 1) -
                                 spline.Dy(i, j - 1))
            );
        }
    });
}

void ReducedAlgorithm::FillDxy(Spline& spline) {
    double threeDivH = 3/xDimension_.H();
    double twelveDivH = 12/xDimension_.H();
    auto unknownsCount = xDimension_.KnotCount();
    auto even = unknownsCount % 2 == 0;
    auto tau = even ? 0 : 2;
    auto eta = even ? -4 : 1;
    auto upsilon = even ? unknownsCount : unknownsCount - 1;
    auto equationsCount = even ? unknownsCount / 2 - 1
                               : unknownsCount / 2;
    int indexes[] = {0, static_cast<int>(spline.ColumnsCount() - 1)};
    for (size_t idx = 0; idx < sizeof(indexes)/sizeof(int); ++idx) {
        auto j = indexes[0];
        auto& tridiagonal = xTridiagonals_.Get();
        auto& rightSideX = tridiagonal.RightSideBuffer();
        for (size_t i = 0; i < equationsCount-1; ++i) {
            auto i21 = 2*(i + 1);
            rightSideX[i] = threeDivH*(spline.Dx(i21 + 2, j) - spline.Dx(i21 - 2, j))
                            - twelveDivH*(spline.Dx(i21 + 1, j) - spline.Dx(i21 - 1, j));
        }
        rightSideX[0] -= spline.Dxy(0, j);
        rightSideX[rightSideX.size() - 1] =
                threeDivH*(spline.Dx(upsilon + tau, j) - spline.Dx(upsilon - 2, j))
                - twelveDivH*(spline.Dx(upsilon + 1, j) - spline.Dx(upsilon - 1, j))
                -eta*spline.Dxy(spline.RowsCount() - 1, j);

        tridiagonal.Solve();

        for (int i = 0; i < rightSideX.size(); ++i) {
            auto i21 = 2*(i + 1);
            spline.SetDx(i21, j, rightSideX[i]);
        }

        for (size_t i = 1; i < spline.RowsCount(); i += 2) {
            spline.SetDxy(i, j,
                         0.25 * (threeDivH * (spline.Dx(i + 1, j) -
                                              spline.Dx(i - 1, j)) -
                                 spline.Dxy(i + 1, j) -
                                 spline.Dxy(i - 1, j))
            );
        }
    }
}

void ReducedAlgorithm::FillDyx(Spline& spline) {
    double threeDivH = 3/yDimension_.H();
    double twelveDivH = 12/yDimension_.H();
    auto unknownsCount = yDimension_.KnotCount();
    auto even = unknownsCount % 2 == 0;
    auto tau = even ? 0 : 2;
    auto eta = even ? -4 : 1;
    auto upsilon = even ? unknownsCount : unknownsCount - 1;
    auto equationsCount = even ? unknownsCount / 2 - 1
                               : unknownsCount / 2;
    utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i) {
//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        auto& tridiagonal = yTridiagonals_.Get();
        auto& rightSideY = tridiagonal.RightSideBuffer();
        for (size_t j = 0; j < equationsCount-1; ++j) {
            auto j21 = 2*(j + 1);
            rightSideY[j] = threeDivH*(spline.Dy(i, j21 + 2) - spline.Dy(i, j21 - 2))
                            - twelveDivH*(spline.Dy(i, j21 + 1) - spline.Dy(i, j21 - 1));
        }
        rightSideY[0] -= spline.Dy(i, 0);
        rightSideY[rightSideY.size() - 1] =
                threeDivH*(spline.Dy(i, upsilon + tau) - spline.Dy(i, upsilon - 2))
                - twelveDivH*(spline.Dy(i, upsilon + 1) - spline.Dy(i, upsilon - 1))
                -eta*spline.Dy(i, spline.ColumnsCount() - 1);

        tridiagonal.Solve();

        for (int j = 0; j < rightSideY.size(); ++j) {
            auto j21 = 2*(j + 1);
            spline.SetDxy(i, j21, rightSideY[j]);
        }

        for (size_t j = 1; j < spline.RowsCount(); j += 2) {
            spline.SetDxy(i, j,
                         0.25 * (threeDivH * (spline.Dy(i, j + 1) -
                                              spline.Dy(i, j - 1)) -
                                 spline.Dxy(i, j + 1) -
                                 spline.Dxy(i, j - 1))
            );
        }
    });
}

void ReducedAlgorithm::Initialize(Spline& spline) {
    InitializeTridiagonals(spline);
    spline.Initialize(function_);
}

bool ReducedAlgorithm::IsParallel() {
    return isParallel_;
}

void ReducedAlgorithm::InParallel(bool value) {
    isParallel_ = value;
}

void ReducedAlgorithm::Parallelize(bool inParallel) {
    xTridiagonals_.Parallelize(inParallel);
    yTridiagonals_.Parallelize(inParallel);
}

void ReducedAlgorithm::InitializeTridiagonals(Spline& spline) {
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().emplace_back(
            Tridiagonal::Factory::CreateReducedTridiagonal(spline.RowsCount()));
    yTridiagonals_.GetAll().emplace_back(
            Tridiagonal::Factory::CreateReducedTridiagonal(spline.ColumnsCount()));
    Parallelize(isParallel_);
}

double ReducedAlgorithm::ExecutionTime() {
    return timer_.ExecutionTime();
}

double ReducedAlgorithm::AllTime() {
    return timer_.AllTime() +
           xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].AllTime();
}
    