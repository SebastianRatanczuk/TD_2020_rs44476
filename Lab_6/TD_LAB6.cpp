//SO IS1 210B LAB06
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
    int A = 2;
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
    int A = 2;
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

    //cout << MaxZa << endl;
    MaxZa -= 3;
    //cout << MaxZa << endl;

    int i;

    for (i = 0; i < len / 2; i++)
    {
        if (array[i] >= MaxZa) break;
    }
    MinID_Za = i;

    //cout << Fk[MinID_Za] << endl;

    for (i = len / 2; i >= 0; i--)
    {
        if (array[i] >= MaxZa) break;
    }
    MaxID_Za = i;

    //cout << Fk[MaxID_Za + 1] << endl;

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

    //cout << "Watek nr " << GetCurrentThreadId() << " Zakonczyl dft" << endl;

    return 0;
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

    N_ = 10;

    for (int i = 0; i < len; i++)
    {
        zA[i] = ASK(x[i], y[i]);
        zF[i] = FSK(x[i], y[i]);
        zP[i] = PSK(x[i], y[i]);
    }

    GenerateData(x, y, len, "Inf");
    GenerateData(x, zA, len, "ZA");
    GenerateData(x, zF, len, "ZF");
    GenerateData(x, zP, len, "ZP");


    HANDLE* threads = new HANDLE[3];
    DWORD* thrdids = new DWORD[3];
    Double_tabe* tables = new Double_tabe[3];


    tables[0].table_Y = zA;
    tables[0].array_length = len;   

    tables[1].table_Y = zF;
    tables[1].array_length = len;    

    tables[2].table_Y = zP;
    tables[2].array_length = len;    

    clock_t t1 = clock();
    for (int i = 0; i < 3; i++)
    {
        threads[i] = CreateThread(NULL, 0, thread, &tables[i], 0, &thrdids[i]);
        if (threads[i] == NULL)
        {
            printf("ERROR");
        }
    }

    for (int i = 0; i < 3; i++)
    {
        long retval;
        WaitForSingleObject(threads[i], INFINITE);        
        CloseHandle(threads[i]);
    }

    clock_t t2 = clock();
    double time = (t2 - t1) / (double)CLOCKS_PER_SEC;
    cout << "Czas trwania dft " << time << endl;

    Complex* Threead_DFT_zA = tables[0].retrived_data;
    Complex* Threead_DFT_zF = tables[1].retrived_data;
    Complex* Threead_DFT_zP = tables[2].retrived_data;    


    double* M_zA1 = new double[len];
    double* M_zF1 = new double[len];
    double* M_zP1 = new double[len];

    double* Mp_zA1 = new double[len];
    double* Mp_zF1 = new double[len];
    double* Mp_zP1 = new double[len];

    double* Fk = new double[len];

    int threshhold = 0;
    for (int i = 0; i < len; i++)
    {
        M_zA1[i] = sqrt(pow(Threead_DFT_zA[i].real(), 2) + pow(Threead_DFT_zA[i].imag(), 2));
        M_zF1[i] = sqrt(pow(Threead_DFT_zF[i].real(), 2) + pow(Threead_DFT_zF[i].imag(), 2));
        M_zP1[i] = sqrt(pow(Threead_DFT_zP[i].real(), 2) + pow(Threead_DFT_zP[i].imag(), 2));

        Mp_zA1[i] = 10 * log10(M_zA1[i]);
        Mp_zF1[i] = 10 * log10(M_zF1[i]);
        Mp_zP1[i] = 10 * log10(M_zP1[i]);

        if (Mp_zA1[i] < threshhold)
            Mp_zA1[i] = 0;
        
        if (Mp_zF1[i] < threshhold)
            Mp_zF1[i] = 0;
        
        if (Mp_zP1[i] < threshhold)
            Mp_zP1[i] = 0;

        Fk[i] = (double)(i)*freq / len;
    }

    GenerateData(Fk, Mp_zA1, len, "Dft_ZA1");
    GenerateData(Fk, Mp_zF1, len, "Dft_ZF1");
    GenerateData(Fk, Mp_zP1, len, "Dft_ZP1");

    double szerokosc_zA = band(Mp_zA1, len, Fk);
    double szerokosc_zF = band(Mp_zF1, len, Fk);
    double szerokosc_zP = band(Mp_zP1, len, Fk);

    cout << "Szerokosc pasma dla zA: " << szerokosc_zA << endl;
    cout << "Szerokosc pasma dla zF: " << szerokosc_zF << endl;
    cout << "Szerokosc pasma dla zP: " << szerokosc_zP << endl;

    system("gnuplot P.plt");
    system("gnuplot P_10.plt");
    system("gnuplot DFT.plt");
}