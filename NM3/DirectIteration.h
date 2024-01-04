#pragma once
#include "Vector.h"
#include "SimpleMatrix.h"

struct trio
{
	double maxValue;
	Vector maxVector;
	int iterationNum;
};

class DirectIteration
{
	int size;
	SimpleMatrix data;

public:
	DirectIteration(int _size);
	DirectIteration(int _size, int border1, int border2, bool print);

	SimpleMatrix getData();
	void printMatrix() const;
	trio findMinEigenValue(double eps, int M);
private:
	double findAngle(Vector& xk, const Vector& xk1);
};