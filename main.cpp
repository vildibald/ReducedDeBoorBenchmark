#include "stdafx.h"
#include "SplineKnots.h"
#include <iostream>
#include <algorithm>
#include "OperationsBenchmark.h"
#include "ComparisonBenchmarkResult.h"
#include <numeric>
#include "StopWatch.h"
#include "FullAlgorithm.h"
#include "ReducedAlgorithm.h"

void MulDivBenchmark() {
    OperationsBenchmark bencher;
    bencher.BenchAll();
}


void PrintSurfaceDeboorResult(ComparisonBenchmarkResult result) {
    std::cout << "Full : " << result.FirstAlg() << std::endl;
    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
    std::cout << "Difference F/R: " << result.Ratio() << std::endl;
}

void SurfaceBenchmark(int numIterations, int numKnots, bool inParallel = false) {
    SurfaceDimension uDim(-20, 20, numKnots);
    SurfaceDimension vDim = uDim;

    MathFunction function = [](double x, double y) {
        return sin(sqrt(x * x + y * y));
    };
    InterpolativeMathFunction f = function;

    FullAlgorithm full(uDim, vDim, f);
    ReducedAlgorithm reduced(uDim, vDim, f);

    std::vector<double> calculatedResults;
    std::vector<double> fullTimes;
    fullTimes.reserve(numIterations);
    std::vector<double> reducedTimes;
    reducedTimes.reserve(numIterations);

    calculatedResults.reserve(numIterations * 2);

    full.InParallel(inParallel);
    reduced.InParallel(inParallel);

    for (size_t i = 0; i < numIterations; i++) {
        auto result = full.Calculate();
        double time = full.ExecutionTime();
        calculatedResults.emplace_back(result.Dxy(1, 1));
        fullTimes.emplace_back(time);
    }

    for (size_t i = 0; i < numIterations; i++) {
        auto result = reduced.Calculate();
        double time = reduced.ExecutionTime();
        calculatedResults.emplace_back(result.Dxy(1, 1));
        reducedTimes.emplace_back(time);
    }

    auto full_time = static_cast<double>(std::accumulate(fullTimes.begin(),
                                                         fullTimes.end(), 0))
                     / static_cast<double>(numIterations);
    auto reduced_time = static_cast<double>(std::accumulate(
            reducedTimes.begin(), reducedTimes.end(), 0))
                        / static_cast<double>(numIterations);
    std::cout << "Ignore " << calculatedResults[0] << std::endl;
    PrintSurfaceDeboorResult(ComparisonBenchmarkResult(full_time, reduced_time));
}


int main() {

    while (true) {
        //std::cout << clock();
        // Console clear ...
        // ... for Windows,
        system("cls");
        // ... for Linux/Unix.
        //system("clear");
        std::cout << "1: Instructions benchmark." << std::endl;
        std::cout << "2: Spline surface benchmark." << std::endl;
        std::cout << "3: Spline surface benchmark (in parallel)." << std::endl;
        std::cout << "Q: End program" << std::endl;
        char input;
        std::cin >> input;
        std::cin.get();
        std::cout << std::endl << "---------------" << std::endl;
        unsigned int num_iterations;
        unsigned int num_knots;
        switch (input) {
            case '1':
                std::cout << "Instructions benchmark" << std::endl << std::endl;
                MulDivBenchmark();
                break;
            case '2':
                std::cout << "Spline surface benchmark" << std::endl << std::endl;
                // At least 10 iterations recomended for stable results.
                std::cout << "Enter number of iterations: " << std::endl;
                std::cin >> num_iterations;
                std::cout << "Enter number of knots: " << std::endl;
                std::cin >> num_knots;
                std::cin.get();
                SurfaceBenchmark(num_iterations, num_knots);
                break;
            case '3':
                std::cout << "Spline surface benchmark (in parallel)" << std::endl << std::endl;
                // At least 10 iterations recomended for stable results.
                std::cout << "Enter number of iterations: " << std::endl;
                std::cin >> num_iterations;
                std::cout << "Enter number of knots: " << std::endl;
                std::cin >> num_knots;
                std::cin.get();
                SurfaceBenchmark(num_iterations, num_knots, true);
                break;

            case 'q':
            case 'Q':
                return 0;
        }

        std::cout << "===================" << std::endl;
        system("pause");
    }
    return 0;

}

