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

	// ��������� ����������� ��������

	Vector eigenValues(size);
	int diapozonValues = border2 - border1 + 1;
	if (diapozonValues < size)
	{
		std::cout << "������������� �������� � ����������� ���������. ������ ��������� ���������: 1,2,.." << size << std::endl;
		for (int i = 0; i < size; i++)
			eigenValues.setValues(i, static_cast<double>(i) + 1);
	}
	else
	{
		for (int i = 0; i < size; i++) // �������� ��������� ����������� ��������
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

	//��������� ����������� vectors

	Vector omega(size, border1, border2); // ��������� �����
	omega = omega / omega.getSecondNorm(); // ������������ �����

	SimpleMatrix E(size, size); // ��������� �������
	for (int i = 0; i < size; i++)
		E.setValue(i, i, 1);

	SimpleMatrix eigenVectors(size, size); // ������� ����������� �������� (������� �����������)
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			eigenVectors.setValue(i, j, E.getValue(i, j) - 2 * omega.getValues(i) * omega.getValues(j));

	SimpleMatrix Lyabda(size, size); // ������������ ������� ������
	for (int i = 0; i < size; i++)
		Lyabda.setValue(i, i, eigenValues.getValues(i));

	data = eigenVectors * Lyabda * eigenVectors;

	// ������ ��� �������� �����������������

	if (print)
	{
		std::cout << "����������� ��������: ";
		eigenValues.print();
		std::cout << "����������� �������:\n";
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

// ����� ���������� �������� �������� ���� ����� ��������� xk � xk1 � ��������� ����

double DirectIteration::findAngle(Vector& xk, const Vector& xk1)
{
	double sum = 0;									//
	for (int i = 0; i < xk.getSize(); i++)			// ��������� ��������� ������������ ��� ����� ������������ ��������������� ���������
		sum += xk.getValues(i) * xk1.getValues(i);	//
	double xkNorm = xk.getSecondNorm();				// 
	xk = xk1; // ��� ����������� ��������

	double answer = sum / (xkNorm * xk1.getSecondNorm()); // ����������� ���� �� ������� 

	return answer;
}

trio DirectIteration::findMinEigenValue(double eps, int M)
{
	int numIteration = 0;		// ����� ����������� ��������
	Vector xk(size, -10, 10);	// ��������� ������ �
	Vector niulast(size, -4, 1);// ������ ����������� ��� ������� ����������� �������
	Vector niuNew(size, 1, 4);	// ������ � �������� ��������� ����������� ������ ������������� �����. ��-�
	double sigma = 0;			// ����� � �������� ��������� ������������ ����������� ��������
	do
	{
		numIteration++;

		niuNew = xk / xk.getSecondNorm();					// ������������ ������� �
		xk = data.solutionHaletskiy(niuNew);				// ��������� �(k+1)

		double sum = 0;
		for (int i = 0; i < size; i++)						// ��������� �����
			sum += niuNew.getValues(i) * xk.getValues(i);	//
		sigma = 1 / sum;									//

	} while (1 - abs(findAngle(niulast, niuNew)) > eps && numIteration < M); // �������� ������� �� ��������� ���� �� ��������� �����. ���� �� ���������� M

	return trio{ sigma, niuNew, numIteration };
}