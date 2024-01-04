#include "DirectIteration.h"
#include <iostream>
#include <iomanip>

DirectIteration::DirectIteration(int _size) : data(_size, _size)
{
	size = _size;
}

DirectIteration::DirectIteration(int _size, int border1, int border2, bool print) : data(_size, _size)
{
	size = _size;

	// получение собственных значений

	Vector eigenValues(size);
	int diapozonValues = border2 - border1 + 1;
	if (diapozonValues < size)
	{
		std::cout << "Повторяющиеся значения в собственных значениях. Выбран фиктивный результат: 1,2,.." << size << std::endl;
		for (int i = 0; i < size; i++)
			eigenValues.setValues(i, static_cast<double>(i) + 1);
	}
	else
	{
		for (int i = 0; i < size; i++) // задаются различные собственные значения
		{
			eigenValues.setValues(i, border1 + static_cast<double>(rand() % diapozonValues));

			auto lyambda = [](const Vector& ptr, int pos, double value)
			{
				for (int i = 0; i < pos; i++)
					if (abs(ptr.getValues(i)) == abs(value)) return true;
				return false;
			};
			while (lyambda(eigenValues, i, eigenValues.getValues(i)))
				eigenValues.setValues(i, border1 + static_cast<double>(rand() % diapozonValues));
		}
	}

	//получение собственных vectors

	Vector omega(size, border1, border2); // получение омеги
	omega = omega / omega.getSecondNorm(); // нормирование омеги

	SimpleMatrix E(size, size); // единичная матрица
	for (int i = 0; i < size; i++)
		E.setValue(i, i, 1);

	SimpleMatrix eigenVectors(size, size); // матрица собственных векторов (матрица Хаусхолдера)
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			eigenVectors.setValue(i, j, E.getValue(i, j) - 2 * omega.getValues(i) * omega.getValues(j));

	SimpleMatrix Lyabda(size, size); // симметричная матрица Лямбда
	for (int i = 0; i < size; i++)
		Lyabda.setValue(i, i, eigenValues.getValues(i));

	data = eigenVectors * Lyabda * eigenVectors;

	// печать для проверки работоспособности

	if (print)
	{
		std::cout << "Собственные значения: ";
		eigenValues.print();
		std::cout << "Собственные вектора:\n";
		eigenVectors.print();
	}
}

SimpleMatrix DirectIteration::getData()
{
	return data;
}

void DirectIteration::printMatrix() const
{
	data.print();
}

// метод возвращает значение косинуса угла между векторами xk и xk1 в радианной мере

double DirectIteration::findAngle(Vector& xk, const Vector& xk1)
{
	double sum = 0;									//
	for (int i = 0; i < xk.getSize(); i++)			// считается скалярное произведение как сумма произведения соответствующих координат
		sum += xk.getValues(i) * xk1.getValues(i);	//
	double xkNorm = xk.getSecondNorm();				// 
	xk = xk1; // для продолжения итерации

	double answer = sum / (xkNorm * xk1.getSecondNorm()); // высчитываем угол по формуле 

	return answer;
}

trio DirectIteration::findMinEigenValue(double eps, int M)
{
	int numIteration = 0;		// число совершённых операций
	Vector xk(size, -10, 10);	// случайный вектор Х
	Vector niulast(size, -4, 1);// вектор необходимый для анализа приближения решения
	Vector niuNew(size, 1, 4);	// вектор к которому стремится собственный вектор максимального собст. зн-я
	double sigma = 0;			// число к которому стремится максимальное собственное значение
	do
	{
		numIteration++;

		niuNew = xk / xk.getSecondNorm();					// нормирование вектора х
		xk = data.solutionHaletskiy(niuNew);				// получение х(k+1)

		double sum = 0;
		for (int i = 0; i < size; i++)						// получение сигмы
			sum += niuNew.getValues(i) * xk.getValues(i);	//
		sigma = 1 / sum;									//

	} while (1 - abs(findAngle(niulast, niuNew)) > eps && numIteration < M); // проверка решение на улучшение либо до заданного епсил. либо до превышения M

	return trio{ sigma, niuNew, numIteration };
}