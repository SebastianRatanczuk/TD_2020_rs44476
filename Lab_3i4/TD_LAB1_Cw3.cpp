#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <time.h>
using namespace std;
typedef complex<double> Complex;

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
    return  1.0 * sin(2. * M_PI * t * 7. + 4. * M_PI);
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

Complex* dft(double*& table, const int& N)
{
    Complex* dfttable = new Complex[N];
    for (int n = 0; n < N; n++)
        for (int k = 0; k < N; k++)
            dfttable[k] += polar(table[n], -2 * M_PI * k * n / N);
    return dfttable;
}

int main(void)
{
    int Tmin = 0;
    int Tmax = 1;
    double Tdelta = 1 / pow(2, 14);


    //calculate X range
    int Interval = Tmax - Tmin;
    //calculate quantity of X instances
    int Quantity = Interval / Tdelta + 1;

    //cout << Quantity << endl; 
    double* tableX = new double[Quantity];
    double* tableY = new double[Quantity];
    double* M = new double[Quantity];
    double* Mp = new double[Quantity];
    double* Fk = new double[Quantity];


    int A = 6, B = 7, C = 4;
    cout << "Preparing done " << Quantity << endl;
    for (int i = 0; i < Quantity; i++)
    {
        tableX[i] = Tdelta * i;
        tableY[i] = S(tableX[i]);
    }
    cout << "Data done" << endl;
    GenerateData(tableX, tableY, Quantity, "daneS.dat");
    cout << "DFT start" << endl;


    int Testsize = Quantity;


    Complex* CDFT = dft(tableY, Testsize);

    cout << "DFT done" << endl;

    for (int i = 0; i < Testsize; i++)
    {
        M[i] = sqrt(pow(CDFT[i].real(), 2) + pow(CDFT[i].imag(), 2));

        Mp[i] = 10 * log10(M[i]);
        if (Mp[i] < 0)
            Mp[i] = 0;
        Fk[i] = (double)i * (1 / Tdelta) / Quantity;
    }
    cout << "Data write Start" << endl;
    GenerateData(tableX, M, Testsize, "daneM.dat");
    GenerateData(tableX, Mp, Testsize, "daneMp.dat");
    GenerateData(tableX, Fk, Testsize, "T3.dat");
    GenerateData(Fk, M, Testsize, "T1.dat");
    GenerateData(Fk, Mp, Testsize, "T2.dat");
    cout << "Data write Done" << endl;

    delete[] tableX;
    delete[] tableY;
    delete[] CDFT;
    return 0;
}