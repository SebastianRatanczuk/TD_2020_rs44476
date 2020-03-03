#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

template<typename T>
void GenerateData(T* Xtable, T* Ytable,int length)
{
    ofstream outfile;
    outfile.open("Data.dat"); 
	outfile << "#\tplot \"Data.dat\" with lines" << endl;

	for (int i = 0; i < length; i++)
	{
		outfile << Xtable[i]<<"\t"<<Ytable[i] << endl;
	}
 
    outfile.close();
}

void Zeros(double A, double B,double C)
{	
	double Delta = pow(B, 2) - 4 * (double)A * C;

	if (Delta < 0)
		cout << "Brak miejsc zerowych" << endl;

	else if (Delta == 0)
	{
		double z1 = B / (-2 * (double)A);
		cout << "Jedno miejsce zerowe " << z1 << endl;
	}
	else
	{
		double z1 = (-B + sqrt(Delta)) / (2 * (double)A);
		double z2 = (-B - sqrt(Delta)) / (2 * (double)A);
		cout << "Dwa miejsce zerowe " << z1 << "\ti " << z2 << endl;
	}
}

int main()
{	
	int Tmin = -10;
	int Tmax = 10;
	double Tdelta = 1. / 100.;

	int Interval = Tmax - Tmin;
	int Quantity = Interval / Tdelta;

	double* tableX = new double[Quantity];
	double* tableY = new double[Quantity];

	unsigned int i = 0;

	//44476
	int A = 6, B = 7, C = 4;

	for (double T = Tmin; T <= Tmax; T += Tdelta)
	{
		tableX[i] = T;
		tableY[i++] = A * pow(T, 2) + B * T + C;
	}

	Zeros(A, B, C);	
	GenerateData(tableX, tableY, Quantity+1);

	return 0;
}