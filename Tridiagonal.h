//
// Created by Viliam on 28.09.2017.
//
#include "stdafx.h"
#include "Timer.h"
#include "MultithreadPreparator.h"
#include <vector>
#include <numeric>
#ifndef SPLINEKNOTS_TRIDIAGONAL_H
#define SPLINEKNOTS_TRIDIAGONAL_H


class Tridiagonal;

class Tridiagonals;

class Tridiagonals final{
    std::vector<Tridiagonal> tridiagonals;

public:
    std::vector<Tridiagonal>& GetAll();

    Tridiagonal& Get();

    void Parallelize(bool inParallel) {
        MultithreadPreparator multithreadPreparator;
        multithreadPreparator.PrepareVector(inParallel, GetAll());
    }

    void Initialize(Tridiagonal tridiagonal);
};

class Tridiagonal final {
    std::vector<double> luBuffer_;
    std::vector<double> rightSideBuffer_;
    double mainDiagonalValue_;
    size_t numEquations_;
    size_t numUnknowns_;


    Tridiagonal(double mainDiagonalValue,
                size_t numEquations,
                size_t numUnknowns
    );

    Timer timer_;

public:


    std::vector<double>& Solve();

    std::vector<double>& ResetBufferAndGet();

    std::vector<double>& Buffer();

    double MainDiagonalValue() const;

    std::vector<double>& RightSideBuffer();

    size_t NumEquations() const;

    size_t NumUnknowns() const;

    double ExecutionTime() {
        return timer_.ExecutionTime();
    }

    double AllTime() {
        return timer_.AllTime();
    }

    class Factory final {
    public:
        static Tridiagonal
        CreateFullTridiagonal(size_t numUnknowns);

        static Tridiagonal
        CreateReducedTridiagonal(size_t numUnknowns);

        static Tridiagonal
        CreateEmptyTridiagonal();

    };

    static double accumulateAllTimes(Tridiagonals& tridiagonals) {
        std::accumulate(tridiagonals.GetAll().begin(), tridiagonals.GetAll().end(),
                        tridiagonals.GetAll()[0].AllTime(),
                        [](double time, Tridiagonal& tridiagonal) {
                            return time + tridiagonal.AllTime();
                        });
    }

    static double accumulateExecutionTimes(Tridiagonals& tridiagonals);
};




#endif //SPLINEKNOTS_TRIDIAGONAL_H
