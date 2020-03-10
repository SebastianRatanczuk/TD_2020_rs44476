#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

template<typename T>
void GenerateData(T* Xtable, T* Ytable, int length, string name)
{
	ofstream outfile;
	outfile.open(name);
	outfile << "#\tplot \"Data.dat\" with lines" << endl;

	for (int i = 0; i < length; i++)
	{
		outfile << Xtable[i] << "\t" << Ytable[i] << endl;
	}

	outfile.close();
}

double S(double t)
{
	return  1.0 * sin(2. * M_PI * t *7. + 4. * M_PI);
}

double Quantization(double x,int depth)
{
	return floor(x *pow(2,depth-1)-1)+pow(2,depth-1)+1;
}

int main()
{
	int Tmin = 0;
	int Tmax = 6;
	double Tdelta = 1. / 1024.;

	//calculate X range
	int Interval = Tmax - Tmin;
	//calculate quantity of X instances
	int Quantity = Interval / Tdelta + 1;

	//cout << Quantity << endl;
	double* tableX = new double[Quantity];
	double* tableY = new double[Quantity];

	//44476
	int A = 6, B = 7, C = 4;

	//prepare X values
	for (int i = 0; i < Quantity; i++)
	{
		tableX[i] = Tdelta * i;
	}

	//generate Y values
	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = S(tableX[i]);
	}
	GenerateData(tableX, tableY, Quantity, "Sin.dat");

	//generate Y values
	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = Quantization(S(tableX[i]),16);
	}
	GenerateData(tableX, tableY, Quantity, "SinQuant.dat");

	delete[] tableX;
	delete[] tableY;

	Tdelta = 1. / 512.;

	Quantity = Interval / Tdelta + 1;

	tableX = new double[Quantity];
	tableY = new double[Quantity];

	for (int i = 0; i < Quantity; i++)
	{
		tableX[i] = Tdelta * i;
	}

	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = Quantization(S(tableX[i]), 8);
	}
	GenerateData(tableX, tableY, Quantity, "SinQuant1.dat");

	return 0;
}