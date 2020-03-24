#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>

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
    int Tmax = 674;
    int Quantity = pow(2, 12);
    int Interval = Tmax - Tmin;
    double Tdelta = (double)Interval / (double)(Quantity - 1);
    
    double* tableX = new double[Quantity];
    double* tableY = new double[Quantity];
    double* M = new double[Quantity];
    
    int A = 6, B = 7, C = 4;

    for (int i = 0; i < Quantity; i++)
    {
        tableX[i] = Tdelta * i;
    }

    for (int i = 0; i < Quantity; i++)
    {
        tableY[i] = S(tableX[i]);
    }

    GenerateData(tableX, tableY,Quantity,"daneS.dat");
    cout << tableY[0] << endl;
    cout << tableY[4095] << endl;       
    
    Complex *CDFT = dft(tableY, Quantity);

    for (int i = 0; i < Quantity; i++)
    {
        cout << CDFT[i] << endl;
    }



    delete[] tableX;
    delete[] tableY;
    delete[] CDFT;
    return 0;
}
