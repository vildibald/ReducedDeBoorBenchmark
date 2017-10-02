#include "stdafx.h"
#include "ComparisonBenchmarkResult.h"


ComparisonBenchmarkResult::ComparisonBenchmarkResult(double firstAlgTime,
                                                     double secondAlgTime)
        : firstAlg_(firstAlgTime), secondAlg_(secondAlgTime){
    ratio_ = firstAlgTime / secondAlgTime;
}

unsigned long long ComparisonBenchmarkResult::FirstAlg() const {
    return firstAlg_;
}

unsigned long long ComparisonBenchmarkResult::SecondAlg() const {
    return secondAlg_;
}

double ComparisonBenchmarkResult::Ratio() const {
    return ratio_;
}
