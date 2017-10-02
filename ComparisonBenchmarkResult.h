#pragma once

class ComparisonBenchmarkResult {
    double firstAlg_;
    double secondAlg_;
    double ratio_;

public:
    ComparisonBenchmarkResult(double firstAlgTime,
                              double secondAlgTime);

    unsigned long long FirstAlg() const;

    unsigned long long SecondAlg() const;

    double Ratio() const;
};
