#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <time.h>
#include <bitset>

using namespace std;

typedef complex<double> Complex;
typedef bitset<8> Bite;

const double kA = 0.4;
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

Bite* convert(string data, bool reverse = false)
{
    Bite* to_ret = new Bite[data.length()];      
    
    for (int i = 0; i < data.length(); i++)
    {
        to_ret[i] = (Bite)(int)data[i];
    }
    
    
      

    return to_ret;
}


double zA(double t,double M_t)
{

    double A1 = 2;
    double A2 = 1;
    double f = 2;

    if (M_t == 1)
    {
        return  A1 * sin(2. * M_PI * f * t);
    }
    else
    {
        return  A2 * sin(2. * M_PI * f * t);
    }
}



int main(void)
{
    double koniec = 2; //w sec
    int czestotliwosc = 8000;

    double Tdelta = 1. / czestotliwosc;
    int ilosc_próbek = koniec * czestotliwosc;

    double* Tablica_Xy = new double[ilosc_próbek];    

    for (int i = 0; i < ilosc_próbek; i++)
    {
        Tablica_Xy[i] = i * Tdelta;
    }   

    string data = "A";
    Bite* bites = convert(data,true);
    cout << *bites << endl;
}