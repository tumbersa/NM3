#pragma once
#include <iostream>
#include "Vector.h"

class SimpleMatrix
{
	int height;
	int width;
	double** data;

public:
	SimpleMatrix();
	SimpleMatrix(int _height, int _width);
	SimpleMatrix(int _height, int _width, int border1, int border2);
	SimpleMatrix(const SimpleMatrix& x);
	~SimpleMatrix();

	SimpleMatrix& operator=(const SimpleMatrix& x);
	SimpleMatrix operator+(const SimpleMatrix& x) const;
	SimpleMatrix operator-(const SimpleMatrix& x) const;
	SimpleMatrix operator*(double x) const;
	Vector operator*(const Vector& x) const;
	SimpleMatrix operator*(const SimpleMatrix& x) const;
	SimpleMatrix transposition() const;
	Vector solutionHaletskiy(const Vector& F) const;

	void setValue(int i, int j, double value);
	double getValue(int i, int j) const;
	void print() const;
};