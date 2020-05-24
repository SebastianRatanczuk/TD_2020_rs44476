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
        binarystring[i] = convert(data[i],true);
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
double ASK(double x, double y,int N,double Td)
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
double FSK(double x, double y,int N, double Td)
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

int main()
{
    double Td = 0.1;
    int Freq = 5000;
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



    /* Generate a new random seed from system time - do this once in your constructor */
    srand(time(0));

    /* Setup constants */
    const static int q = 15;
    const static float c1 = (1 << q) - 1;
    const static float c2 = ((int)(c1 / 3)) + 1;
    const static float c3 = 1.f / c1;

    /* random number in range 0 - 1 not including 1 */
    float random = 0.f;

    /* the white noise */
    float noise = 0.f;

    for (int i = 0; i < len; i++)
    {
        random = ((float)rand() / (float)(RAND_MAX + 1));
        noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * (c3*1);        
        zA[i] += noise;
        zF[i] += noise;
        zP[i] += noise;
    }

    GenerateData(x, y, len, "Inf");
    GenerateData(x, zA, len, "ZA");
    GenerateData(x, zF, len, "ZF");
    GenerateData(x, zP, len, "ZP");


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
    double H = 0.2;
    for (int i = 0; i < len; i++)
    {
        if (P[i] < H)
            MA[i] = 0;
        else
            MA[i] = 1;
    }

    GenerateData(x, P, len, "Calka_zA");
    GenerateData(x, MA, len, "Sygnal_zA");


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

    int packets = recived_fA.length() / 7;

    string* reshapePackets = new string [packets];


    //TUTAJ wybór 
    for (int i = 0; i < recived_fA.length(); i++)
    {       
        reshapePackets[i / 7] += recived_fP[i];
    }


    cout << endl << "Odzyskane pakiety" << endl;
    for (int i = 0; i < packets; i++)
    {
        cout << reshapePackets[i] << endl;
    }    

    cout << endl << "Zdekodowane pakiety" << endl;
    string* DecodedHamming = new string[packets];
    for (int i = 0; i < packets; i++)
    {
        int * tmp = DecodeHamming(reshapePackets[i]);
        string tmpstr;

        for (int j = 0; j < 7; j++)
        {
            tmpstr += to_string(tmp[j]);
        }

        DecodedHamming[i] = tmpstr;
        cout << DecodedHamming[i] << endl;
    }
    
    cout << endl << "Wydobycie danych" << endl;
    string* PacketsHamming = new string[packets];

    for (int i = 0; i < packets; i++)
    {
        string tmp;
        tmp += DecodedHamming[i][2];
        tmp += DecodedHamming[i][4];
        tmp += DecodedHamming[i][5];
        tmp += DecodedHamming[i][6];
        PacketsHamming[i]=tmp;
        cout << PacketsHamming[i] << endl;
    }

    cout << endl << "Laczenie pakietow danych" << endl;
    string* HammingData = new string[packets/2];
    for (int i = 0; i < packets/2; i++)
    {
        HammingData[i] += PacketsHamming[2 * i];
        HammingData[i] += PacketsHamming[2 * i + 1];

        cout << HammingData[i] << endl;
    }

    cout << endl << "Rozszyfrowanie wiadomosci" << endl;

    for (int i = 0; i < packets / 2; i++)
    {
        cout << char(Binto10(HammingData[i]));
    }
}