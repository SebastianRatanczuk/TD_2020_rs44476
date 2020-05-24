//sebastian ratanczuk
//sebastian-ratanczuk@zut.edu.pl

#define _USE_MATH_DEFINES

#include <windows.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>
#include <time.h>
#include <bitset>

using namespace std;

typedef complex<double> Complex;
typedef bitset<8> Bite;

int N_ = 1;
double Td_ = 1;

struct Data
{
    int length;
    int bite;
    double* x;
    double* y;
};

void GenerateData(Data data, string name)
{
    ofstream outfile;
    outfile.open(name);
    outfile << "#\tplot \"Data.dat\" with lines" << endl;

    for (int i = 0; i < data.length; i++)
    {
        outfile << data.x[i] << "\t" << data.y[i] << endl;
    }

    outfile.close();
}

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

Complex* dft(double* tab, const int N) {

    Complex* dtf = new Complex[N];
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < N; k++) {
            double m = -2 * M_PI * k * n / N;
            dtf[k] += polar(tab[n], m);
        }
    }
    return dtf;
}

void reverseStr(string& str)
{
    int n = str.length();

    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}

string convert(char data, bool reverse = false)
{
    Bite bite = (Bite)(int)data;
    string reverseorder = bite.to_string();
    if (reverse)
    {
        reverseStr(reverseorder);
        return reverseorder;
    }
    return bite.to_string();
}

string convert(string data, bool reverse = false)
{
    string to_ret = "";

    for (int i = 0; i < data.length(); i++)
    {
        to_ret += convert(data[i], reverse);
    }

    return to_ret;
}

Data generateSignal(string bite, int Sample_freq)
{
    int bite_length = Td_ * Sample_freq;

    int word_length = bite_length * bite.length();

    double* x = new double[word_length];
    double* signal = new double[word_length];

    int currentpos = 0;
    for (int i = 0; i < bite.length(); i++)
    {
        int bit = bite[i] - '0';

        for (currentpos; currentpos < (i + 1) * bite_length; currentpos++)
        {
            signal[currentpos] = bit;
            x[currentpos] = currentpos * 1. / Sample_freq;

            if (currentpos > word_length)
            {
                cout << "ERROR 1";
            }
        }
    }

    Data data;
    data.length = word_length;
    data.x = x;
    data.y = signal;
    data.bite = bite_length;

    return data;
}

double ASK(double x, double y)
{
    int A1 = 3;
    int A2 = 1;
    int f = N_ / Td_;
    if (y == 1)
        return  A1 * sin(2. * M_PI * f * x);
    else if (y == 0)
        return  A2 * sin(2. * M_PI * f * x);
    else
        return  -101;
}

double FSK(double x, double y)
{
    int A = 3;
    int f1 = (N_ * 10) / Td_;
    int f2 = (N_) / Td_;
    if (y == 1)
        return  A * sin(2. * M_PI * f1 * x);
    else if (y == 0)
        return  A * sin(2. * M_PI * f2 * x);
    else
        return  -101;
}

double PSK(double x, double y)
{
    int A = 3;
    int f = N_ / Td_;
    if (y == 1)
        return  A * sin(2. * M_PI * f * x);
    else if (y == 0)
        return  A * sin(2. * M_PI * f * x + M_PI);
    else
        return  -101;
}


double band(double* array, int len, double* Fk)
{
    double MaxZa = 0;
    int MinID_Za = 0, MaxID_Za = 0;
    for (int i = 0; i < len / 2; i++)
    {
        if (MaxZa < array[i]) MaxZa = array[i];
    }

    MaxZa -= 3;

    int i;

    for (i = 0; i < len / 2; i++)
    {
        if (array[i] >= MaxZa) break;
    }
    MinID_Za = i;

    for (i = len / 2; i >= 0; i--)
    {
        if (array[i] >= MaxZa) break;
    }
    MaxID_Za = i;

    return Fk[MaxID_Za + 1] - Fk[MinID_Za];
}

struct Double_tabe {

    double* table_Y;
    int array_length;
    Complex* retrived_data;
};

DWORD WINAPI thread(LPVOID data)
{
    Double_tabe* dane = (Double_tabe*)data;

    double* tab = dane->table_Y;
    const int N = dane->array_length;

    Complex* dtf = new Complex[N];
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < N; k++) {
            double m = -2 * M_PI * k * n / N;
            dtf[k] += polar(tab[n], m);
        }
    }

    dane->retrived_data = dtf;

    return 0;
}

double calka(double* funkcja, int poczatek, int koniec, double h)
{
    double suma = 0;
    for (int i = poczatek; i < koniec; i++)
    {
        suma += funkcja[i];
    }

    return h * suma;
}

int main(void)
{
    string test = "Tajne";
    string bites = convert(test, false);
    int dlogoscnapisu = test.length();

    Td_ = 0.1;
    int freq = 8000;

    Data data = generateSignal(bites, freq);

    double* x = data.x;
    double* y = data.y;
    int len = data.length;
    int bitelen = data.bite;

    double* zA = new double[len];
    double* zF = new double[len];
    double* zP = new double[len];

    N_ = 2;

    for (int i = 0; i < len; i++)
    {
        zA[i] = ASK(x[i], y[i]);
        zF[i] = FSK(x[i], y[i]);
        zP[i] = PSK(x[i], y[i]);
    }

    GenerateData(x, y, data.bite * 10, "Inf_10");
    GenerateData(x, zA, data.bite * 10, "ZA_10");
    GenerateData(x, zF, data.bite * 10, "ZF_10");
    GenerateData(x, zP, data.bite * 10, "ZP_10");
   

    GenerateData(x, y, len, "Inf");
    GenerateData(x, zA, len, "ZA");
    GenerateData(x, zF, len, "ZF");
    GenerateData(x, zP, len, "ZP");

    
    double* X = new double[len];
    double* X2 = new double[len];
    double* P = new double[len];
    double* M = new double[len];

    //zA
    int f = N_ / Td_;
    //tworzenie sygnalu X
    for (int i = 0; i < len; i++)
    {
        X[i]=zA[i]* 3 * sin(2. * M_PI * f * x[i]);
    }

    GenerateData(x, X, len, "XzA");

    //calkowanie    
    for (int i = 0; i < dlogoscnapisu*8; i++)
    {
        double n = calka(X, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        for (int j = bitelen * i; j < bitelen * (i + 1); j++)
        {
            P[j] = n;
        }
    }

    // cyfryzacja?
    double H = 0.2;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            M[i] = 0;
        else
            M[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zA");
    GenerateData(x, M, len, "Sygnal_zA");


    //zP

    for (int i = 0; i < len; i++)
    {
        X[i] = zP[i] * 3 * sin(2. * M_PI * f * x[i]);
    }

    GenerateData(x, X, len, "XzP");

    //calkowanie    
    for (int i = 0; i < dlogoscnapisu * 8; i++)
    {
        double n = calka(X, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        for (int j = bitelen * i; j < bitelen * (i + 1); j++)
        {
            P[j] = n;
        }
    }

    // cyfryzacja?
    H = 0;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            M[i] = 0;
        else
            M[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zP");
    GenerateData(x, M, len, "Sygnal_zP");

    //zF
    
    int f1 = (N_ * 10) / Td_;
    int f2 = (N_) / Td_;

    //tworzenie sygnalu X1 i X2
    for (int i = 0; i < len; i++)
    {
        X[i] = zF[i] * 3 * sin(2. * M_PI * f1 * x[i]);
        X2[i] = zF[i] * 3 * sin(2. * M_PI * f2* x[i]);
    }

    GenerateData(x, X, len, "XzF");
    GenerateData(x, X2, len, "XzF2");

    //calkowanie X 
    for (int i = 0; i < dlogoscnapisu * 8; i++)
    {
        double n = calka(X, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        double n2 = calka(X2, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        //cout << n << " " << n2 << endl;
        for (int j = bitelen * i; j < bitelen * (i + 1); j++)
        {
            P[j] = n-n2;
        }
    }

    // cyfryzacja?
    H = 0;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            M[i] = 0;
        else
            M[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zF");
    GenerateData(x, M, len, "Sygnal_zF");
}