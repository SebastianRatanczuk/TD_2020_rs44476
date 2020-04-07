#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <time.h>
using namespace std;
typedef complex<double> Complex;

const double kA = 0.5;
const double kP = 1;


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

Complex* dft(double*& table, const int& N, const int& K)
{
    Complex* dfttable = new Complex[K];
    Complex i(0.0, 1);

    Complex Wn = exp(i * 2. * M_PI / (double)N);

    for (int k = 0; k < K; k++)
    {
        dfttable[k] = 0;
        for (int n = 0; n < N; n++)
        {
            dfttable[k] += table[n] * pow(Wn, -k * n);
        }
    }

    return dfttable;
}

double* idft(Complex*& X, const int& N)
{
    double* double_table = new double[N];
    Complex i(0.0, 1);
    Complex Wn = exp(i * 2. * M_PI / (double)N);
    for (int n = 0; n < N; n++)
    {
        double_table[n] = 0;
        for (int k = 0; k < N; k++)
        {
            Complex tmp = pow(Wn, k * n) * X[k];
            double_table[n] += tmp.real();
        }
        double_table[n] *= 1. / N;
    }
    return double_table;
}

double Fm(double t, int N)
{
    double Sum = 0;

    for (int n = 1; n <= N; n++)
    {
        Sum += (cos(12. * t * pow(n, 2)) + cos(16. * t * n)) / pow(n, 2);
    }

    return Sum;
}

double Fm(double t)
{
    return  1. * sin(2. * M_PI * t * 4. + 0. * M_PI);
}

double Fn(double t)
{   
    return  1. * sin(2. * M_PI * t * 800.*2. + 0. * M_PI);
}

double M(double t)
{
    return  2. * sin(2. * M_PI * Fm(t,2));
}

double zA(double t)
{
    return  (kA * M(t) + 1.) * cos(2. * M_PI * Fn(t));
}

double zP(double t)
{
    return  cos(2. * M_PI*Fn(t)+kP*M(t));
}

int main(void)
{
    int koniec = 1; //w sec
    int czestotliwosc = 1000000;

    double Tdelta = 1. / czestotliwosc;
    int ilosc_próbek = 1 * czestotliwosc;

    double*  Tablica_Xy = new double[ilosc_próbek];
    double* Tablica_Za = new double[ilosc_próbek];
    double* Tablica_Zp = new double[ilosc_próbek];
    double* Tablica_M = new double[ilosc_próbek];
    double* Tablica_S = new double[ilosc_próbek];

    for (int i = 0; i < ilosc_próbek; i++)
    {
        Tablica_Xy[i] = i * Tdelta;
        Tablica_Za[i] = zA(Tablica_Xy[i]);
        Tablica_Zp[i] = zP(Tablica_Xy[i]);
        Tablica_M[i] = M(Tablica_Xy[i]);
        Tablica_S[i] = Fm(Tablica_Xy[i],2);
    }

    GenerateData(Tablica_Xy, Tablica_Za, ilosc_próbek,"Za.dat");
    GenerateData(Tablica_Xy, Tablica_Zp, ilosc_próbek,"Zp.dat");
    GenerateData(Tablica_Xy, Tablica_M, ilosc_próbek,"M.dat");
    GenerateData(Tablica_Xy, Tablica_S, ilosc_próbek,"S.dat");

    return 0;
}