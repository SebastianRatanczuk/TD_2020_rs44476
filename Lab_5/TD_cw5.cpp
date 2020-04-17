#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <time.h>
using namespace std;
typedef complex<double> Complex;

const double kA = 0.8;
const double kP = 0.9;

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

double M(double t)
{
    return  1. * sin(2. * M_PI * Fm(t,2) * t);
}

double zA(double t)
{
    return  (kA * M(t) + 1.) * cos(2. * M_PI * 800. * t);
    //return  cos(2. * M_PI * 800. * t);
}

double zP(double t)
{
    return  cos(2. * M_PI * 800 * t + kP * M(t));
}

int main(void)
{
    double koniec = 1; //w sec
    int czestotliwosc = 8000;

    double Tdelta = 1. / czestotliwosc;
    int ilosc_próbek = koniec * czestotliwosc;

    double* Tablica_Xy = new double[ilosc_próbek];
    double* Tablica_Za = new double[ilosc_próbek];
    double* Tablica_Zp = new double[ilosc_próbek];    
    double* M_Za = new double[ilosc_próbek];
    double* Mp_Za = new double[ilosc_próbek];
    double* M_Zp = new double[ilosc_próbek];
    double* Mp_Zp = new double[ilosc_próbek];
    double* Fk = new double[ilosc_próbek];

    for (int i = 0; i < ilosc_próbek; i++)
    {
        Tablica_Xy[i] = i * Tdelta;
        Tablica_Za[i] = zA(Tablica_Xy[i]);
        Tablica_Zp[i] = zP(Tablica_Xy[i]);       
    }

    GenerateData(Tablica_Xy, Tablica_Za, ilosc_próbek,"Za.dat");
    GenerateData(Tablica_Xy, Tablica_Zp, ilosc_próbek,"Zp.dat");
     

    int ilosc_harmo = 8000;
    Complex* Tdft_Za = dft(Tablica_Za, ilosc_próbek, ilosc_harmo);
    Complex* Tdft_Zp = dft(Tablica_Zp, ilosc_próbek, ilosc_harmo);

    for (int i = 0; i < ilosc_harmo; i++)
    {
        M_Za[i] = sqrt(pow(Tdft_Za[i].real(), 2) + pow(Tdft_Za[i].imag(), 2));
        M_Zp[i] = sqrt(pow(Tdft_Zp[i].real(), 2) + pow(Tdft_Zp[i].imag(), 2));
        
       
        Mp_Za[i] = 10 * log10(M_Za[i]);
        if (Mp_Za[i] < 0)
            Mp_Za[i] = 0;
        
        Mp_Zp[i] = 10 * log10(M_Zp[i]);
        if (Mp_Zp[i] < 0)
            Mp_Zp[i] = 0;  

        Fk[i] = (double)(i)*czestotliwosc / ilosc_próbek;
    }


    GenerateData(Fk, M_Za, ilosc_harmo, "MZA.dat");
    GenerateData(Fk, Mp_Za, ilosc_harmo, "MPZA.dat");
    GenerateData(Fk, M_Zp, ilosc_harmo, "MZP.dat");
    GenerateData(Fk, Mp_Zp, ilosc_harmo, "MPZP.dat");

    return 0;
}