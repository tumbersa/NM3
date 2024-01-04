#include <iostream>
#include <iomanip>
#include "SimpleMatrix.h"

SimpleMatrix::SimpleMatrix()
{
	height = 0;
	width = 0;
	data = nullptr;
}

SimpleMatrix::SimpleMatrix(int _height, int _width)
{
	height = _height;
	width = _width;

	data = new double* [height];
	for (int i = 0; i < height; i++)
		data[i] = new double[width];

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			data[i][j] = 0;
}

SimpleMatrix::SimpleMatrix(int _height, int _width, int border1, int border2)
{
	height = _height;
	width = _width;

	data = new double* [height];
	for (int i = 0; i < height; i++)
		data[i] = new double[width];

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			data[i][j] = border1 + static_cast<double>(rand() % (border2 - border1 + 1));
}

SimpleMatrix::SimpleMatrix(const SimpleMatrix& x)
{
	height = x.height;
	width = x.width;

	data = new double* [height];
	for (int i = 0; i < height; i++)
		data[i] = new double[width];

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			data[i][j] = x.data[i][j];
}

SimpleMatrix::~SimpleMatrix()
{
	for (int i = 0; i < height; i++)
		delete[] data[i];
	delete[] data;
}

SimpleMatrix& SimpleMatrix::operator=(const SimpleMatrix& x)
{
	if (&x != this)
	{
		height = x.height;
		width = x.width;

		data = new double* [height];
		for (int i = 0; i < height; i++)
			data[i] = new double[width];

		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++)
				data[i][j] = x.data[i][j];
	}
	return *this;
}

SimpleMatrix SimpleMatrix::operator+(const SimpleMatrix& x) const
{
	if (height != x.height || width != x.width)
	{
		std::cout << "Разные размерности. Сложение матриц невозможно.\n";
		return *this;
	}
	SimpleMatrix answ(height, width);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			answ.data[i][j] = data[i][j] + x.data[i][j];

	return answ;
}

SimpleMatrix SimpleMatrix::operator-(const SimpleMatrix& x) const
{
	if (height != x.height || width != x.width)
	{
		std::cout << "Разные размерности. Вычитание матриц невозможно.\n";
		return *this;
	}
	SimpleMatrix answ(height, width);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			answ.data[i][j] = data[i][j] - x.data[i][j];

	return answ;
}

SimpleMatrix SimpleMatrix::operator*(double x) const
{
	SimpleMatrix answ(height, width);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			answ.data[i][j] = data[i][j] * x;

	return answ;
}

Vector SimpleMatrix::operator*(const Vector& x) const
{
	if (width != x.getSize())
	{
		std::cout << "Разные размерности. Умножение матрицы на вектор невозможно.\n";
		return x;
	}
	Vector answ(height);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			answ.setValues(i, answ.getValues(i) + data[i][j] * x.getValues(j));

	return answ;
}

SimpleMatrix SimpleMatrix::operator*(const SimpleMatrix& x) const
{
	if (width != x.height)
	{
		std::cout << "Разные размерности. Умножение матриц невозможно.\n";
		return *this;
	}
	SimpleMatrix answ(height, x.width);

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < x.width; j++)
		{
			double sum = 0;
			for (int k = 0; k < width; k++)
				sum += data[i][k] * x.data[k][j];
			answ.data[i][j] = sum;
		}
	}
	return answ;
}

SimpleMatrix SimpleMatrix::transposition() const
{
	if (width != height)
	{
		std::cout << "Матрица не квадратная. транспонирование невозможно.\n";
		return *this;
	}
	SimpleMatrix answ(height, width);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			answ.data[i][j] = data[j][i];

	return answ;
}

Vector SimpleMatrix::solutionHaletskiy(const Vector& F) const
{
	if (width != height)
	{
		std::cerr << "Ширина и длина не совпадают. Решения системы отсутствуют\n";
		return Vector(height);
	}

	//создание матрицы BC
	double** BC = new double* [height];
	for (int i = 0; i < height; i++)
		BC[i] = new double[width];

	//получение матрицы BC
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < height; j++)
		{
			double sum = 0;
			if (j <= i)//считаем B
			{
				for (int k = 0; k < j; k++)
					sum += BC[i][k] * BC[k][j];

				BC[i][j] = data[i][j] - sum;
			}
			else //считаем C
			{
				for (int k = 0; k < i; k++)
					sum += BC[i][k] * BC[k][j];

				BC[i][j] = (data[i][j] - sum) / BC[i][i];
			}
		}
	}

	Vector X(height);
	Vector Y(height);

	// Подсчёт Y
	for (int i = 0; i < height; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += BC[i][j] * Y.getValues(j);

		Y.setValues(i, (F.getValues(i) - sum) / BC[i][i]);
	}

	// Подсчёт X
	for (int i = height - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < height; j++)
			sum += BC[i][j] * X.getValues(j);

		X.setValues(i, Y.getValues(i) - sum);
	}

	//удаление BC
	for (int i = 0; i < height; i++)
		delete[] BC[i];
	delete[] BC;

	return X;
}

void SimpleMatrix::setValue(int i, int j, double value)
{
	data[i][j] = value;
}

double SimpleMatrix::getValue(int i, int j) const
{
	return data[i][j];
}

void SimpleMatrix::print() const
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
			std::cout << std::setw(12) << std::setprecision(5) << data[i][j];
		std::cout << std::endl << std::endl;
	}
}