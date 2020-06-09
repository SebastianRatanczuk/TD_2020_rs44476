#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>
#include <time.h>
#include <bitset>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#include <errno.h>





using namespace std;

typedef complex<double> Complex;
typedef bitset<8> Bite;



struct HammingData
{
    string data;
    string* binarystring;
    string* packets;
    string* Hammingpackets;

    int PacketsAmmount;
};


//Struk0tura przechowojaca wartosci funkcji
struct Data
{
    int length;
    int bite;
    double* x;
    double* y;
};


//Klasa Macierzy
class Matrix
{
public:
    int** matrix;
    int x;
    int y;

    Matrix()
    {
        matrix = NULL;
        x = 0;
        y = 0;
    }

    Matrix(int x, int y)
    {
        this->x = x;
        this->y = y;

        matrix = new int* [x];
        for (int i = 0; i < x; i++)
            matrix[i] = new int[y];
    }

    Matrix(int x, int y, string data)
    {
        this->x = x;
        this->y = y;

        if (x * y != data.length())
            throw ("Data length error");

        matrix = new int* [x];
        for (int i = 0; i < x; i++)
            matrix[i] = new int[y];


        int count = 0;
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
                matrix[i][j] = data[count++] - '0';
        }
    }

    Matrix(int x, int y, int* data, int len)
    {
        this->x = x;
        this->y = y;

        if (x * y != len)
            throw ("Data length error");

        matrix = new int* [x];
        for (int i = 0; i < x; i++)
            matrix[i] = new int[y];


        int count = 0;
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
                matrix[i][j] = data[count++];
        }
    }


    //Wypisuje Macierz na ekran
    void show()
    {
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
                cout << matrix[i][j];
        }
        cout << endl;
    }


    //Mnozenie Macierzowe
    Matrix operator*(Matrix b)
    {
        if (this->y == b.x)
        {
            Matrix to_ret(this->x, b.y);

            for (int i = 0; i < to_ret.x; i++)
            {
                for (int j = 0; j < to_ret.y; j++)
                {
                    int suma = 0;

                    for (int h = 0; h < b.x; h++)
                        suma = suma + this->matrix[i][h] * b.matrix[h][j];

                    to_ret.matrix[i][j] = suma;
                }
            }
            return to_ret;
        }
        else throw("Dimmensions dont match");
    }


    //Modulo dla kazdego elementu macierzy 
    void operator%(int t)
    {
        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++)
                matrix[i][j] = matrix[i][j] % t;
    }


    //Rozwiniecie maciezy do jednowymiarowej tablicy
    int* list()
    {
        int* to_ret = new int[x * y];
        int c = 0;
        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++)
                to_ret[c++] = matrix[i][j];
        return to_ret;
    }


    //Rozszezenie tablicy bez wpisywania danych
    void append(int dane)
    {
        this->x++;
        int** new_matrix = new int* [x];
    }
};


//Wypisanie danych do pliku
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


//Wypisanie danych do pliku
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


//Negacja wskazanego bitu
void NEG(int* data, int number)
{
    data[number] = !data[number];
}


//Funkcja odwracajaca string
void reverseStr(string& str)
{
    int n = str.length();

    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}


//Zmiana string na ascii binarnie 
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


//Zmiana Z sys binarnego do dziesietnego
int Binto10(int* data, int range)
{
    int to_ret = 0;

    for (int i = 0; i < range; i++)
        to_ret += data[i] * pow(2, i);


    return to_ret;
}


//Zmiana z sys binarnego do dziesietnego
int Binto10(string data)
{
    int to_ret = 0;

    for (int i = 0; i < data.length(); i++)
        to_ret += (data[i] - '0') * pow(2, i);

    return to_ret;
}


//Kodowanie hamminga na pojednczym pakiecie
int* CodeHamming(string Data)
{
    string gstr = "1101101110000111010000100001";
    Matrix G(7, 4, gstr);

    Matrix D(4, 1, Data);

    Matrix K = G * D;
    K % 2;

    return K.list();
}


//Dekodwanie Haminga dla pojedynczego pakietu
int* DecodeHamming(string Data)
{
    string hstr = "101010101100110001111";
    Matrix H(3, 7, hstr);

    Matrix Kp(7, 1, Data);

    Matrix S = H * Kp;
    S % 2;
    //S.show();

    //cout << "Pozycja blednego bitu " << endl;
    int x = Binto10(S.list(), 3);
    //cout << x << endl;

    int* to_ret = Kp.list();

    if (x != 0)
    {
        cout << "Blad w pakiecie" << endl;
        NEG(to_ret, x - 1);

    }

    //cout << endl << "Naprawiony Kod" << endl;
    /*for (int i = 0; i < 7; i++)
        cout << to_ret[i];
    cout << endl;*/

    return to_ret;
}


//Kodowanie Hamminga dla podanego napisu
HammingData EncodeHamming(string data)
{
    string* binarystring = new string[data.length()];

    cout << "Znaki zamienione na asci" << endl;

    for (int i = 0; i < data.length(); i++)
    {
        binarystring[i] = convert(data[i], true);
        cout << binarystring[i] << endl;
    }



    string* packets = new string[data.length() * 2];

    for (int i = 0; i < data.length(); i++)
    {
        string tmp = binarystring[i];
        packets[2 * i] = tmp.substr(0, tmp.length() / 2);
        packets[2 * i + 1] = tmp.substr(tmp.length() / 2);
    }



    cout << endl << "Bajty podzielone na pakiety" << endl;

    for (int i = 0; i < data.length() * 2; i++)
    {
        cout << packets[i] << endl;
    }



    string* HammingPackets = new string[data.length() * 2];

    cout << endl << "Pakiety zabezpieczone kodowaniem Hamminga" << endl;

    for (int i = 0; i < data.length() * 2; i++)
    {
        int* tmp = CodeHamming(packets[i]);
        string tmpstr;

        for (int j = 0; j < 7; j++)
            tmpstr += to_string(tmp[j]);

        HammingPackets[i] = tmpstr;

        cout << HammingPackets[i] << endl;
    }
    cout << endl;



    HammingData to_ret;
    to_ret.data = data;
    to_ret.binarystring = binarystring;
    to_ret.packets = packets;
    to_ret.Hammingpackets = HammingPackets;
    to_ret.PacketsAmmount = data.length() * 2;

    return to_ret;
}


//Rozwiniecie do pojedynczego stringa
string Unroll(string* data, int size)
{
    string to_ret;

    for (int i = 0; i < size; i++)
        to_ret += data[i];

    return to_ret;
}


//Generowanie ciagu binarnego
Data generateSignal(string bite, int Sample_freq, double Td)
{
    int bite_length = Td * Sample_freq;

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


//Kluczowanie ASK
double ASK(double x, double y, int N, double Td)
{
    int A1 = 3;
    int A2 = 1;
    int f = N / Td;
    if (y == 1)
        return  A1 * sin(2. * M_PI * f * x);
    else if (y == 0)
        return  A2 * sin(2. * M_PI * f * x);
    else
        return  -101;
}


//Kluczowanie FSK
double FSK(double x, double y, int N, double Td)
{
    int A = 3;
    int f1 = (N * 10) / Td;
    int f2 = (N) / Td;
    if (y == 1)
        return  A * sin(2. * M_PI * f1 * x);
    else if (y == 0)
        return  A * sin(2. * M_PI * f2 * x);
    else
        return  -101;
}


//Kluczowanie PSK
double PSK(double x, double y, int N, double Td)
{
    int A = 3;
    int f = N / Td;
    if (y == 1)
        return  A * sin(2. * M_PI * f * x);
    else if (y == 0)
        return  A * sin(2. * M_PI * f * x + M_PI);
    else
        return  -101;
}


//Calkowanie numeryczne metoda trapezow
double calka(double* funkcja, int poczatek, int koniec, double h)
{
    double suma = 0;
    for (int i = poczatek; i < koniec; i++)
    {
        suma += funkcja[i];
    }

    return h * suma;
}


//Wszystkie tablice funkcji kluczjacych
struct Signals
{
    int len;

    double* zA;
    double* zF;
    double* zP;
    double* x;
};


//Funkcja kluczjaca po wszystkich metodach
Signals SendSignal(Data Signal, double Td, int N)
{
    int len = Signal.length;
    double* x = Signal.x;
    double* y = Signal.y;

    double* zA = new double[len];
    double* zF = new double[len];
    double* zP = new double[len];


    for (int i = 0; i < len; i++)
    {
        zA[i] = ASK(x[i], y[i], N, Td);
        zF[i] = FSK(x[i], y[i], N, Td);
        zP[i] = PSK(x[i], y[i], N, Td);
    }

    Signals to_ret;

    to_ret.len = len;
    to_ret.x = x;
    to_ret.zA = zA;
    to_ret.zF = zF;
    to_ret.zP = zP;

    return to_ret;
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


void dec(string recived)
{
    

    int packets = recived.length() / 7;
    string* reshapePackets = new string[packets];
    //TUTAJ wybór 
    for (int i = 0; i < recived.length(); i++)
    {
        reshapePackets[i / 7] += recived[i]; //TU ZMIEN
    }


    /*cout << endl << "Odzyskane pakiety" << endl;
    for (int i = 0; i < packets; i++)
    {
        cout << reshapePackets[i] << endl;
    }*/

    //cout << endl << "Zdekodowane pakiety" << endl;
    string* DecodedHamming = new string[packets];
    for (int i = 0; i < packets; i++)
    {
        int* tmp = DecodeHamming(reshapePackets[i]);
        string tmpstr;

        for (int j = 0; j < 7; j++)
        {
            tmpstr += to_string(tmp[j]);
        }

        DecodedHamming[i] = tmpstr;
        //cout << DecodedHamming[i] << endl;
    }

   //cout << endl << "Wydobycie danych" << endl;
    string* PacketsHamming = new string[packets];

    for (int i = 0; i < packets; i++)
    {
        string tmp;
        tmp += DecodedHamming[i][2];
        tmp += DecodedHamming[i][4];
        tmp += DecodedHamming[i][5];
        tmp += DecodedHamming[i][6];
        PacketsHamming[i] = tmp;
        //cout << PacketsHamming[i] << endl;
    }

    //cout << endl << "Laczenie pakietow danych" << endl;
    string* HammingData = new string[packets / 2];
    for (int i = 0; i < packets / 2; i++)
    {
        HammingData[i] += PacketsHamming[2 * i];
        HammingData[i] += PacketsHamming[2 * i + 1];

        //cout << HammingData[i] << endl;
    }    

    for (int i = 0; i < packets / 2; i++)
    {
        cout << char(Binto10(HammingData[i]));
    }

    cout << endl;
}

int main()
{
    double Td = 0.1;
    int Freq = 1000;
    string data = "Ala ma kota";
    cout << "Zdanie do wyslania" << endl << data << endl;

    HammingData Hamm = EncodeHamming(data);
    string bites = Unroll(Hamm.Hammingpackets, Hamm.PacketsAmmount);
    cout << "Bity przygotowane do nadania w jedym ciagu" << endl;
    cout << bites << endl;

    Data Signal = generateSignal(bites, Freq, Td);

    int N = 2;
    Signals Recived = SendSignal(Signal, Td, N);

    double* x = Recived.x;
    double* y = Signal.y;

    double* zA = Recived.zA;
    double* zF = Recived.zF;
    double* zP = Recived.zP;

    int bitelen = Signal.bite;

    int len = Signal.length;


    srand(time(0));    
    const static int q = 15;
    const static float c1 = (1 << q) - 1;
    const static float c2 = ((int)(c1 / 3)) + 1;
    const static float c3 = 1.f / c1;    
    float random = 0.f;    
    float noise = 0.f;
    

    cout << "Wypis1" << endl;
    GenerateData(x, y, len, "Inf");
    GenerateData(x, zA, len, "ZA");
    GenerateData(x, zF, len, "ZF");
    GenerateData(x, zP, len, "ZP");
    
    double* zA_Noise = new double[len];
    double* zF_Noise = new double[len];
    double* zP_Noise = new double[len];

    double alf = 0.2;
    for (int i = 0; i < len; i++)
    {
        random = ((float)rand() / (float)(RAND_MAX + 1));
        noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * (c3 * 1);

        zA_Noise[i] = (zA[i] * alf) + (noise * (1. - alf));
        zF_Noise[i] = (zF[i] * alf) + (noise * (1. - alf));
        zP_Noise[i] = (zP[i] * alf) + (noise * (1. - alf));
    }
      
    
    cout << "Wypis2" << endl;
   
    GenerateData(x, zA_Noise, len, "ZA_Noise");
    GenerateData(x, zF_Noise, len, "ZF_Noise");
    GenerateData(x, zP_Noise, len, "ZP_Noise");



    /*cout << "Dft start" << endl;

    HANDLE* threads = new HANDLE[6];
    DWORD* thrdids = new DWORD[6];
    Double_tabe* tables = new Double_tabe[6];


    tables[0].table_Y = zA;
    tables[0].array_length = len;

    tables[1].table_Y = zF;
    tables[1].array_length = len;

    tables[2].table_Y = zP;
    tables[2].array_length = len;

    tables[3].table_Y = zA_Noise;
    tables[3].array_length = len;

    tables[4].table_Y = zF_Noise;
    tables[4].array_length = len;

    tables[5].table_Y = zP_Noise;
    tables[5].array_length = len;


    clock_t t1 = clock();
    for (int i = 0; i < 6; i++)
    {
        threads[i] = CreateThread(NULL, 0, thread, &tables[i], 0, &thrdids[i]);
        if (threads[i] == NULL)
        {
            printf("ERROR");
        }
    }

    for (int i = 0; i < 6; i++)
    {
        long retval;
        WaitForSingleObject(threads[i], INFINITE);
        CloseHandle(threads[i]);
    }

    clock_t t2 = clock();
    double time = (t2 - t1) / (double)CLOCKS_PER_SEC;
    cout << "Czas trwania dft " << time << endl;

    Complex* Thread_DFT_zA = tables[0].retrived_data;
    Complex* Thread_DFT_zF = tables[1].retrived_data;
    Complex* Thread_DFT_zP = tables[2].retrived_data;

    Complex* Thread_DFT_zA_N = tables[3].retrived_data;
    Complex* Thread_DFT_zF_N = tables[4].retrived_data;
    Complex* Thread_DFT_zP_N = tables[5].retrived_data;


    double* M_zA1 = new double[len];
    double* M_zF1 = new double[len];
    double* M_zP1 = new double[len];

    double* Mp_zA1 = new double[len];
    double* Mp_zF1 = new double[len];
    double* Mp_zP1 = new double[len];

    double* M_zA1N = new double[len];
    double* M_zF1N = new double[len];
    double* M_zP1N = new double[len];

    double* Mp_zA1N = new double[len];
    double* Mp_zF1N = new double[len];
    double* Mp_zP1N = new double[len];

    double* Fk = new double[len];

    int threshhold = 0;
    for (int i = 0; i < len; i++)
    {
        M_zA1[i] = sqrt(pow(Thread_DFT_zA[i].real(), 2) + pow(Thread_DFT_zA[i].imag(), 2));
        M_zF1[i] = sqrt(pow(Thread_DFT_zF[i].real(), 2) + pow(Thread_DFT_zF[i].imag(), 2));
        M_zP1[i] = sqrt(pow(Thread_DFT_zP[i].real(), 2) + pow(Thread_DFT_zP[i].imag(), 2));

        Mp_zA1[i] = 10 * log10(M_zA1[i]);
        Mp_zF1[i] = 10 * log10(M_zF1[i]);
        Mp_zP1[i] = 10 * log10(M_zP1[i]);

        if (Mp_zA1[i] < threshhold)
            Mp_zA1[i] = 0;

        if (Mp_zF1[i] < threshhold)
            Mp_zF1[i] = 0;

        if (Mp_zP1[i] < threshhold)
            Mp_zP1[i] = 0;

        M_zA1N[i] = sqrt(pow(Thread_DFT_zA_N[i].real(), 2) + pow(Thread_DFT_zA_N[i].imag(), 2));
        M_zF1N[i] = sqrt(pow(Thread_DFT_zF_N[i].real(), 2) + pow(Thread_DFT_zF_N[i].imag(), 2));
        M_zP1N[i] = sqrt(pow(Thread_DFT_zP_N[i].real(), 2) + pow(Thread_DFT_zP_N[i].imag(), 2));

        Mp_zA1N[i] = 10 * log10(M_zA1N[i]);
        Mp_zF1N[i] = 10 * log10(M_zF1N[i]);
        Mp_zP1N[i] = 10 * log10(M_zP1N[i]);

        if (Mp_zA1N[i] < threshhold)
            Mp_zA1N[i] = 0;

        if (Mp_zF1N[i] < threshhold)
            Mp_zF1N[i] = 0;

        if (Mp_zP1N[i] < threshhold)
            Mp_zP1N[i] = 0;

        Fk[i] = (double)(i)*Freq / len;
    }

    GenerateData(Fk, Mp_zA1, len, "Dft_ZA");
    GenerateData(Fk, Mp_zF1, len, "Dft_ZF");
    GenerateData(Fk, Mp_zP1, len, "Dft_ZP");

    GenerateData(Fk, Mp_zA1N, len, "Dft_ZAN");
    GenerateData(Fk, Mp_zF1N, len, "Dft_ZFN");
    GenerateData(Fk, Mp_zP1N, len, "Dft_ZPN");*/

    zA = zA_Noise;
    zF = zF_Noise;
    zP = zP_Noise;

    int dlogoscnapisu = Hamm.PacketsAmmount * 7;

    double* X = new double[len];
    double* X2 = new double[len];
    double* P = new double[len];
    double* MA = new double[len];
    double* MF = new double[len];
    double* MP = new double[len];

    //zA    
    int f = N / Td;
    //tworzenie sygnalu X
    for (int i = 0; i < len; i++)
    {
        X[i] = zA[i] * 3 * sin(2. * M_PI * f * x[i]);
    }

    GenerateData(x, X, len, "XzA");

    //calkowanie    
    for (int i = 0; i < dlogoscnapisu; i++)
    {
        double n = calka(X, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        for (int j = bitelen * i; j < bitelen * (i + 1); j++)
        {
            P[j] = n;
        }
    }

    // cyfryzacja?
    double H = 0.1;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            MA[i] = 0;
        else
            MA[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zA");
    GenerateData(x, MA, len, "Sygnal_zA");
    cout << "Za done" << endl;


    //zP

    for (int i = 0; i < len; i++)
    {
        X[i] = zP[i] * 3 * sin(2. * M_PI * f * x[i]);
    }

    GenerateData(x, X, len, "XzP");

    //calkowanie    
    for (int i = 0; i < dlogoscnapisu; i++)
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
            MP[i] = 0;
        else
            MP[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zP");
    GenerateData(x, MP, len, "Sygnal_zP");
    cout << "Zp done" << endl;

    //zF

    int f1 = (N * 10) / Td;
    int f2 = (N) / Td;

    //tworzenie sygnalu X1 i X2
    for (int i = 0; i < len; i++)
    {
        X[i] = zF[i] * 3 * sin(2. * M_PI * f1 * x[i]);
        X2[i] = zF[i] * 3 * sin(2. * M_PI * f2 * x[i]);
    }

    GenerateData(x, X, len, "XzF");
    GenerateData(x, X2, len, "XzF2");

    //calkowanie X 
    for (int i = 0; i < dlogoscnapisu; i++)
    {
        double n = calka(X, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        double n2 = calka(X2, bitelen * i, bitelen * (i + 1), x[1] - x[0]);
        //cout << n << " " << n2 << endl;
        for (int j = bitelen * i; j < bitelen * (i + 1); j++)
        {
            P[j] = n - n2;
        }
    }

    // cyfryzacja?
    H = 0;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            MF[i] = 0;
        else
            MF[i] = 1;
    }

    cout << "Zf done" << endl;

    GenerateData(x, P, len, "Calka_zF");
    GenerateData(x, MF, len, "Sygnal_zF");

    string recived_fA;
    string recived_fF;
    string recived_fP;
    for (int i = 0; i < Hamm.PacketsAmmount * 7; i++)
    {
        recived_fA += MA[bitelen * i] + '0';
        recived_fF += MF[bitelen * i] + '0';
        recived_fP += MP[bitelen * i] + '0';
    }

    cout << "Odebrane dane kluczowania ASK" << endl << recived_fA << endl;
    cout << "Odebrane dane kluczowania FSK" << endl << recived_fF << endl;
    cout << "Odebrane dane kluczowania PSK" << endl << recived_fP << endl;    
        
    int zaber = 0;
    int zfber = 0;
    int zpber = 0;
    for (int i = 0; i < recived_fA.length(); i++)
    {
         
        if (bites[i] != recived_fA[i])
            zaber++;
        if (bites[i] != recived_fF[i])
            zfber++;
        if (bites[i] != recived_fP[i])
            zpber++;
    }

    
    cout << endl << "Ask ";
    dec(recived_fA);
    cout << zaber << " " << (double)zaber / bites.length() << endl;

    cout << "Fsk ";
    dec(recived_fF);
    cout << zfber << " " << (double)zfber / bites.length() << endl;

    cout << "Psk ";
    dec(recived_fP);
    cout << zpber << " " << (double)zpber / bites.length() << endl;    
}