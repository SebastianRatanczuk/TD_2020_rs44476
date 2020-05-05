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

int main()
{
    int Freq = 1000;
    double TDCLK = 0.05;
    int CLK_AM = 33;
    int CLK_LEN = CLK_AM * TDCLK * Freq;
    int Bit_len = TDCLK * Freq;
    double* X = new double[CLK_LEN];

    double* CLK = new double[CLK_LEN];
    double* TTL = new double[CLK_LEN];
    double* MANC = new double[CLK_LEN];
    double* NRZI = new double[CLK_LEN];
    double* BAMI = new double[CLK_LEN];

    int flag = 1;

    //clk

    for (int i = 0; i < CLK_AM; i++)
    {
        for (int j = i * Bit_len; j < (i + 1) * Bit_len; j++)
        {
            CLK[j] = flag;
            X[j] = (double)j / Freq;
        }
        flag = !flag;
    }

    GenerateData(X, CLK, CLK_LEN, "CLK");

    string bit = "1010110010001100";


    //ttl
    flag = 0;
    int a = 0;
    for (int i = 1; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] < CLK[i])
        {
            if (a < bit.length())
                flag = bit[a++] - '0';
            else
                flag = 0;
        }
        TTL[i] = flag;
    }

    GenerateData(X, TTL, CLK_LEN, "TTL");

    //manchaster
    flag = 0;
    for (int i = 1; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] < CLK[i])//rosnacy
        {            
            if (TTL[i - 1] == TTL[i])
            {
                if (TTL[i] == 1)
                    flag = -1;
                else
                    flag = 1;
            }
        }

        if (CLK[i - 1] > CLK[i])//malejacy
        {
            if (TTL[i] == 1)
                flag = 1;
            else
                flag = -1;
        }

        MANC[i] = flag;
    }

    GenerateData(X, MANC, CLK_LEN, "Manchester");

    //nrzi
    flag = 0;
    for (int i = 1; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] > CLK[i])//malejacy
        {
            if (TTL[i] == 1)
            {
                if (flag == 0)
                    flag = 1;
                else
                    flag = flag * -1;
            }
                
        }

        NRZI[i] = flag;
    }

    GenerateData(X, NRZI, CLK_LEN, "NRZI");

    //bami

    flag = 0;
    int f1 = 1;
    for (int i = 1; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] < CLK[i])//rosnacy
        {
            if (TTL[i] == 1)
            {
                flag = f1;
                f1 = f1 * -1; ;
            }
            else
                flag = 0;
        }
        BAMI[i] = flag;
    }

    GenerateData(X, BAMI, CLK_LEN, "BAMI");

}