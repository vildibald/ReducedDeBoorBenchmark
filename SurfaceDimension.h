#pragma once

struct SurfaceDimension {
    double min, max;
    int knotCount;
    double h;


    SurfaceDimension(double min, double max, int knotCount);

    double H() const;

    double Min() const;

    double Max() const;

    int KnotCount() const;

};

