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
		outfile << Xtable[i]<<"\t"<<Ytable[i] << endl;
	}
 
    outfile.close();
}

double X(double t)
{
	return  6. * pow(t, 2) + 7. * t + 4.;	
}

double Y(double t)
{
	return 2. * pow(X(t), 2) + 12. * cos(t);
}

double Z(double t)
{
	return sin(2. * M_PI * 7. * t) * X(t) - 0.2 * log10(abs(Y(t) + M_PI));
}

double U(double t)
{
	return sqrt(abs(pow(Y(t), 2) * Z(t))) - 1.8 * sin(0.4 * t * Z(t) * X(t));
}

double V(double t)
{
	if (t >= 0)
	{
		if (1 >= t)
		{
			if (t >= 0.7)
			{
				return pow(t, -0.662) + 0.77 * sin(8 * t);
			}
			else if (0.22 <= t)
			{
				return 0.63 * t * sin(125 * t);
			}
			else if (0.22 > t)
			{
				return (1 - 7 * t) * sin((2 * M_PI + t + 10) / (t + 0.4));
			}
			else
				return -1002;
		}
		else 
			return -1001;
	}
	else
		return -1000;
}

double P(double t, int N)
{
	double Sum = 0;

	for (int n = 1; n <= N; n++)
	{
		Sum += (cos(12 * t * pow(n, 2)) + cos(16 * t * n)) / pow(n, 2);
	}

	return Sum;
}

int main()
{	
	int Tmin = 0;
	int Tmax = 1;
	double Tdelta = 1. / 22050.;

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
		tableX[i] = Tdelta*i;		
	}

	//generate Y values
	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = Y(tableX[i]);
	}
	GenerateData(tableX, tableY, Quantity,"Y.dat");


	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = Z(tableX[i]);
	}
	GenerateData(tableX, tableY, Quantity, "Z.dat");


	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = U(tableX[i]);
	}
	GenerateData(tableX, tableY, Quantity, "U.dat");


	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = V(tableX[i]);
	}
	GenerateData(tableX, tableY, Quantity, "V.dat");
	
	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = P(tableX[i],2);
	}
	GenerateData(tableX, tableY, Quantity, "P2.dat");
		

	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = P(tableX[i],4);
	}
	GenerateData(tableX, tableY, Quantity, "P4.dat");


	for (int i = 0; i < Quantity; i++)
	{
		tableY[i] = P(tableX[i], A * 10 + B);
	}
	GenerateData(tableX, tableY, Quantity, "PAB.dat");


	return 0;
}