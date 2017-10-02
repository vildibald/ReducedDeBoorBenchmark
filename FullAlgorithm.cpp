#include "FullAlgorithm.h"
#include "StopWatch.h"
#include "utils.h"

FullAlgorithm::FullAlgorithm(const SurfaceDimension xDimension,
                             const SurfaceDimension yDimension,
                             const InterpolativeMathFunction f)
        : xDimension_(xDimension), yDimension_(yDimension), xTridiagonals_(), yTridiagonals_(),
          function_(f), isParallel_(false) {
}

Spline FullAlgorithm::Calculate() {
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

void FullAlgorithm::FillDx(Spline& spline) {
    double threeDivH = 3 / xDimension_.H();
    utils::For(0, static_cast<int>(spline.ColumnsCount()), 1, isParallel_, [&](int j){
//    for (size_t j = 0; j < spline.ColumnsCount(); ++j) {
        auto& tridiagonal = xTridiagonals_.Get();
        auto& rightSideX = tridiagonal.RightSideBuffer();
        for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
            rightSideX[i - 1] = threeDivH * (spline.Z(i + 1, j) - spline.Z(i - 1, j));
        }
        rightSideX[0] -= spline.Dx(0, j);
        rightSideX[rightSideX.size() - 1] -= spline.Dx(spline.RowsCount() - 1, j);
        tridiagonal.Solve();
        for (int i = 0; i < rightSideX.size(); ++i) {
            spline.SetDx(i + 1, j, rightSideX[i]);
        }
    });
}

void FullAlgorithm::FillDy(Spline& spline) {
    double threeDivH = 3 / yDimension_.H();
    utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i){
//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        auto& tridiagonal = yTridiagonals_.Get();
        auto& rightSideY = tridiagonal.RightSideBuffer();
        for (size_t j = 1; j < spline.RowsCount() - 1; ++j) {
            rightSideY[j - 1] = threeDivH * (spline.Z(i, j + 1) - spline.Z(i, j - 1));
        }
        rightSideY[0] -= spline.Dy(i, 0);
        rightSideY[rightSideY.size() - 1] -= spline.Dy(i, spline.ColumnsCount() - 1);
        tridiagonal.Solve();
        for (int j = 0; j < rightSideY.size(); ++j) {
            spline.SetDy(i, j + 1, rightSideY[j]);
        }
    });
}

void FullAlgorithm::FillDxy(Spline& spline) {
    auto& tridiagonal = xTridiagonals_.Get();
    auto& rightSideX = tridiagonal.RightSideBuffer();
    double threeDivH = 3 / xDimension_.H();
    int j = 0;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideX[i - 1] = threeDivH * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideX[0] -= spline.Dxy(0, j);
    rightSideX[rightSideX.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    tridiagonal.Solve();
    for (int i = 0; i < rightSideX.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideX[i]);
    }

    j = spline.ColumnsCount() - 1;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideX[i - 1] = threeDivH * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideX[0] -= spline.Dxy(0, j);
    rightSideX[rightSideX.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    tridiagonal.Solve();
    for (int i = 0; i < rightSideX.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideX[i]);
    }
}

void FullAlgorithm::FillDyx(Spline& spline) {
    double threeDivH = 3 / yDimension_.H();
    utils::For(0, static_cast<int>(spline.RowsCount()), 1, isParallel_, [&](int i){
//    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        auto& tridiagonal = yTridiagonals_.Get();
        auto& rightSideY = tridiagonal.RightSideBuffer();
        for (size_t j = 1; j < spline.RowsCount() - 1; ++j) {
            rightSideY[j - 1] = threeDivH * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1));
        }
        rightSideY[0] -= spline.Dxy(i, 0);
        rightSideY[rightSideY.size() - 1] -= spline.Dxy(i, spline.ColumnsCount() - 1);
        tridiagonal.Solve();
        for (int j = 0; j < rightSideY.size(); ++j) {
            spline.SetDxy(i, j + 1, rightSideY[j]);
        }
    });
}

void FullAlgorithm::Initialize(Spline& spline) {
    InitializeTridiagonals(spline);
    spline.Initialize(function_);
}

bool FullAlgorithm::IsParallel() {
    return isParallel_;
}

void FullAlgorithm::InParallel(bool value) {
    isParallel_ = value;
}

void FullAlgorithm::Parallelize(bool inParallel) {
    xTridiagonals_.Parallelize(inParallel);
    yTridiagonals_.Parallelize(inParallel);
}

void FullAlgorithm::InitializeTridiagonals(Spline& spline) {
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().emplace_back(
            Tridiagonal::Factory::CreateFullTridiagonal(spline.RowsCount()));
    yTridiagonals_.GetAll().emplace_back(
            Tridiagonal::Factory::CreateFullTridiagonal(spline.ColumnsCount()));
    Parallelize(isParallel_);
}

double FullAlgorithm::ExecutionTime() {
    return timer_.ExecutionTime();
}

double FullAlgorithm::AllTime() {
    return timer_.AllTime() +
           xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].AllTime();
}
    