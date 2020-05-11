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
    int CLK_AM = 34;
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
        if (CLK[i - 1] < CLK[i]) //rosnacy
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


    cout <<"Sygnal pierwotny:\t"<< bit << endl;

    //dekodowanie ttl 
    double* D_TTL = new double[bit.length()];
    a = 0;
    for (int i = 1; i < CLK_LEN; i++) 
    {
        if (CLK[i - 1] < CLK[i])//rosnacy
        {
            D_TTL[a++] = TTL[i];
        }
    }


    cout << "TTL:\t\t\t";
    for (int i = 0; i < bit.length(); i++)
    {
        cout << D_TTL[i];
    }
    cout << endl;

    //dekodowanie manshester 
    double* D_MAN = new double[bit.length()];
    a = 0;
    double* NEW_CLK = new double[CLK_LEN];
    int offset = Bit_len * 0.25;
    for (int i = 0; i < CLK_LEN; i++)
    {
        NEW_CLK[i] = (i-offset < 0) ? 1 : CLK[i-offset]; 
    }

    GenerateData(X, NEW_CLK, CLK_LEN, "OFFSET");

    for (int i = Bit_len*2; i < CLK_LEN; i++) 
    {
        if (NEW_CLK[i - 1] > NEW_CLK[i])//malejacy
        {
            D_MAN[a++] = (MANC[i] < 0) ? 0 : 1;
            //D_MAN[a++] = MANC[i];
        }
    }
    cout << "Manchester:\t\t";

    for (int i = 0; i < bit.length(); i++)
    {
        cout << D_MAN[i];
    }
    cout << endl;

    //dekodowanie nrzi
    double* D_NRZItmp = new double[bit.length()]; 
    double* D_NRZI = new double[bit.length()]; 
    a = 0;

    for (int i = 1; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] > CLK[i])//male          
        {
            D_NRZItmp[a++] = (NRZI[i] < 0) ? 0 : 1;
        }
    }

    cout << "NRZI:\t\t\t";

    for (int i = 0; i < bit.length(); i++)
    {
        D_NRZI[i] = (int) D_NRZItmp[i] ^ (int)D_NRZItmp[i+1];
    }

    for (int i = 0; i < bit.length(); i++)
    {
        cout << D_NRZI[i];
    }
    cout << endl;


    //bami
    double* D_BAMI = new double[bit.length()];
    a = 0;    

    for (int i = 2; i < CLK_LEN; i++)
    {
        if (CLK[i - 1] < CLK[i])         
        {
            if (BAMI[i] == 1)
                D_BAMI[a++] = 1;
            if (BAMI[i] == -1)
                D_BAMI[a++] = 1;
            if (BAMI[i] == 0)
                D_BAMI[a++] = 0;
        }
    }

    cout << "BAMI:\t\t\t";
    for (int i = 0; i < bit.length(); i++)
    {
        cout << abs(D_BAMI[i]);
    }
    cout << endl;

}