#include "stdafx.h"
#include "SurfaceDimension.h"

SurfaceDimension::SurfaceDimension(double min, double max, int knotCount)
        : min(min),
          max(max),
          knotCount(knotCount),
          h(abs(max - min)/(knotCount - 1)) {
}

double SurfaceDimension::H() const {
    return h;
}

double SurfaceDimension::Min() const {
    return min;
}

double SurfaceDimension::Max() const {
    return max;
}

int SurfaceDimension::KnotCount() const {
    return knotCount;
}

