#include "stdafx.h"
#include "Spline.h"
#include "utils.h"
#include <iostream>

void Spline::Print()
{
	using namespace std;
	cout << "---------- Knot matrix ----------" << endl;
	for (size_t i = 0; i < xDimension_.knotCount; i++)
	{
		cout << "Row " << i << " :\n";
		for (size_t j = 0; j < yDimension_.knotCount; j++)
		{
			cout << j << ":\n"
				<< "z: " << z[i][j] << '\n'
				<< "dx: " << dx[j][i] << '\n'
				<< "dy: " << dy[i][j] << '\n'
				<< "dxy: " << dxy[i][j] << '\n';
		}
		cout << endl;
	}
	cout << "-------------------------------" << endl;
}
//
Spline::Spline()
	:xDimension_(0,0,0),yDimension_(0,0,0),
    z(), dx(), dy(), dxy()
{
}

Spline Spline::Null()
{
	Spline nullval;
	return nullval;
}

bool Spline::IsNull()
{
    if (xDimension_.knotCount < 1 || yDimension_.knotCount < 1)
		return true;
	return false;
}

Spline::Spline(SurfaceDimension rowDimension, SurfaceDimension columnDimension)
	: xDimension_(rowDimension), yDimension_(columnDimension)
, z(), dx(), dy(), dxy()
{
	z.resize(xDimension_.knotCount);
	dx.resize(xDimension_.knotCount);
	dy.resize(xDimension_.knotCount);
	dxy.resize(xDimension_.knotCount);

	for (int i = 0; i < xDimension_.knotCount; ++i) {
		z[i].resize(yDimension_.knotCount);
		dx[i].resize(yDimension_.knotCount);
		dy[i].resize(yDimension_.knotCount);
		dxy[i].resize(yDimension_.knotCount);
	}
}

void Spline::Initialize(InterpolativeMathFunction mathFunction) {
	auto hx = xDimension_.H();
	auto hy = yDimension_.H();
	auto u = xDimension_.min;
	// Init Z
	for (auto i = 0; i < xDimension_.knotCount; i++, u += hx) {
		auto v = xDimension_.min;
		for (auto j = 0; j < yDimension_.knotCount; j++, v += hy) {
			auto z = mathFunction.Z()(u, v);
			SetZ(i, j, z);
		}
	}
	// Init Dx
	u = xDimension_.min;
	auto uKnotCountMin1 = xDimension_.knotCount - 1;
	for (auto j = 0; j < yDimension_.knotCount; j++) {
//            auto v = vDimension.min;
		auto x0 = X(0);
		auto y0 = Y(0);
		auto x1 = X(uKnotCountMin1);
		auto y1 = Y(uKnotCountMin1);
		auto d0 = mathFunction.Dx()(x0, y0);
		auto d1 = mathFunction.Dx()(x1, y1);
		SetDx(0, j, d0);
		SetDx(uKnotCountMin1, j, d1);
	}
	// Init Dy
	auto vKnotCountMin1 = yDimension_.knotCount - 1;
	for (auto i = 0; i < xDimension_.knotCount; i++) {
		auto x0 = X(0);
		auto y0 = Y(0);
		auto x1 = X(vKnotCountMin1);
		auto y1 = Y(vKnotCountMin1);
		auto d0 = mathFunction.Dy()(x0, y0);
		auto d1 = mathFunction.Dy()(x1, y1);
		SetDy(i, 0, d0);
		SetDy(i, vKnotCountMin1, d1);

	}
	// Init Dxy
	SetDxy(0, 0, mathFunction.Dxy()(X(0), Y(0)));
	SetDxy(uKnotCountMin1, 0, mathFunction.Dxy()(X(0),
														Y(uKnotCountMin1)));
	SetDxy(0, vKnotCountMin1, mathFunction.Dxy()(X(vKnotCountMin1),
														Y(0)));
	SetDxy(uKnotCountMin1, vKnotCountMin1, mathFunction.Dxy()(X(vKnotCountMin1),
																	 Y(uKnotCountMin1)));
}
