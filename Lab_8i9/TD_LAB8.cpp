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

int main()
{

    string test = "DU";
    string bites = convert(test, false);
    string werbit = "1010110010001100";
    //bites = werbit;
    int* x = new int[34];
    int* clock = new int[34];
    int* ttl = new int[34];
    int* man = new int[34];
    int* nr = new int[34];

    for (int i = 0; i < 34; i++)
    {
        x[i] = i;
        clock[i] = (i+1) % 2;    
        ttl[i] = 0;
        man[i] = 0;
        nr[i] = 0;
    }

    for (int i = 0; i < 16; i++)
    {
        ttl[2*i] = bites[i]-'0';
        ttl[2*i + 1] = bites[i] - '0';
    }

    for (int i = 0; i < 16; i++)
    {
        if (bites[i] - '0' == 0)
        {
            man[2 * i] = 1;
            man[2 * i + 1] = -1;
        }
        else if (bites[i] - '0' == 1)
        {
            man[2 * i] = -1;
            man[2 * i + 1] = 1;
        }
        else
        {
            man[2 * i] = -2;
            man[2 * i + 1] = -2;
        }
    }

    GenerateData(x, clock, 34, "Clk");
    GenerateData(x, ttl, 34, "Ttl");
    GenerateData(x, man, 34, "Man");
    
    cout << bites << endl;
}
