#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>


std::vector<double>&
Tridiagonal::ResetBufferAndGet() {
    auto& buffer = luBuffer_;
    std::fill(buffer.begin(), buffer.end(), 1);
    return buffer;
}

std::vector<double>&
Tridiagonal::Buffer() {
    return luBuffer_;
}

std::vector<double>&
Tridiagonal::Solve() {
    auto& buffer = Buffer();
    utils::solveDeboorTridiagonalSystemBuffered(mainDiagonalValue_,
                                                &rightSideBuffer_.front(), numEquations_,
                                                &buffer.front());
    return rightSideBuffer_;
}

Tridiagonal::Tridiagonal(double lhs1Coeficient,
                         size_t numEquations,
                         size_t numUnknowns)
        : mainDiagonalValue_(lhs1Coeficient),
          luBuffer_(numEquations),
          rightSideBuffer_(numEquations),
          numEquations_(numEquations),
          numUnknowns_(numUnknowns),
          timer_(Timer()) {
    luBuffer_.assign(numEquations, 0);
    rightSideBuffer_.assign(numEquations, 0);
}


std::vector<double>& Tridiagonal::RightSideBuffer() {
    return rightSideBuffer_;
}

size_t Tridiagonal::NumEquations() const {
    return numEquations_;
}

size_t Tridiagonal::NumUnknowns() const {
    return numUnknowns_;
}


double Tridiagonal::MainDiagonalValue() const {
    return mainDiagonalValue_;
}

Tridiagonal
Tridiagonal::Factory::CreateFullTridiagonal(size_t numUnknowns) {
    Tridiagonal tridiagonal(
            4,
            numUnknowns - 2,
            numUnknowns
    );
    return tridiagonal;
}

Tridiagonal
Tridiagonal::Factory::CreateReducedTridiagonal(size_t numUnknowns) {
    auto even = numUnknowns%2 == 0;
    auto numEquations = even ? numUnknowns/2 - 2 : numUnknowns/2 - 1;
    Tridiagonal tridiagonal(
            -14,
            numEquations,
            numUnknowns
    );
    return tridiagonal;
}

Tridiagonal Tridiagonal::Factory::CreateEmptyTridiagonal() {
    Tridiagonal tridiagonal(
            0,
            0,
            0
    );
    return tridiagonal;
}

std::vector<Tridiagonal>& Tridiagonals::GetAll() {
    return tridiagonals;
}

Tridiagonal& Tridiagonals::Get() {
    if(tridiagonals.size()==1){
        return tridiagonals[0];
    }
    return tridiagonals[omp_get_thread_num()];
}

void Tridiagonals::Initialize(Tridiagonal tridiagonal) {
    GetAll().clear();
    GetAll().emplace_back(std::move(tridiagonal));
}
