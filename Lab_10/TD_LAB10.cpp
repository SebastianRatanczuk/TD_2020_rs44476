#include <iostream>
#include <string>

using namespace std;

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
            matrix[i] = new int [y];
        

        int count = 0;
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)            
                matrix[i][j] = data[count++] - '0';                        
        }
    }

    Matrix(int x, int y, int* data,int len)
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

    void show()
    {
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)            
                cout << matrix[i][j]<< " ";
            
            cout << endl;
        }
        cout << endl;
    }

    Matrix operator*(Matrix b)
    {
        if (this->y == b.x)
        {
            Matrix to_ret(this->x,b.y);

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

    void operator%(int t)
    {
        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++)
                matrix[i][j] = matrix[i][j]%t;
    }

    int* list()
    {
        int* to_ret = new int[x * y];
        int c = 0;
        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++)
                to_ret[c++] = matrix[i][j];
        return to_ret;
    }

    void append(int dane)
    {
        this->x++;
        int** new_matrix = new int* [x];
    }
};


void NEG(int* data,int number)
{
    data[number] = !data[number];
}

int Binto10(int* data, int range)
{
    int to_ret = 0;

    for (int i = 0; i < range; i++)    
        to_ret += data[i] * pow(2, i);
    

    return to_ret;
}


int Parry(int* data, int len)
{
    int to_ret = 0;
    for (int i = 0; i < len; i++)    
        to_ret += data[i];
    
    return to_ret % 2;
}

int* CodeHamming(string Data)
{
    string gstr = "1101101110000111010000100001";
    Matrix G(7, 4, gstr);

    Matrix D(4, 1, Data);

    Matrix K = G * D;
    K % 2;
    K.show();

    return K.list();
}

int* DecodeHamming(int* Data)
{
    string hstr = "101010101100110001111";
    Matrix H(3, 7, hstr);

    Matrix Kp(7, 1, Data, 7);

    Matrix S = H * Kp;
    S % 2;
    S.show();

    cout <<"Pozycja blednego bitu "<< endl;
    int x = Binto10(S.list(), 3);
    cout << x << endl;

    int* to_ret = Kp.list();

    if(x!=0)
        NEG(to_ret, x - 1);

    cout << endl << "Naprawiony Kod" << endl;
    for (int i = 0; i < 7; i++)
        cout << to_ret[i] << endl;

    return to_ret;
}


int* CodeSECDED(string data)
{
    cout << "Kod hamminga dla lewej czesci" << endl;
    int* k = CodeHamming(data);
    int* SEC = new int[8];

    for (int i = 0; i < 7; i++)    
       SEC[i] = k[i];
    
    SEC[7] = Parry(k,7);
    cout << "Kod SECDED dla lewej czesci" << endl;
    for (int i = 0; i < 8; i++)
        cout << SEC[i] << endl;
    return SEC;
}

int* DecodeSECDED(int *h)
{   
    
    int p4 = Parry(h, 7);    
    if (p4 != h[7])    
        cout << "Wykryto blad proba naprawy" << endl;    
    
    int p1 = (h[0] + h[2] + h[4] + h[6]) % 2;
    int p2 = (h[1] + h[2] + h[5] + h[6]) % 2;
    int p3 = (h[3] + h[4] + h[5] + h[6]) % 2;

    int n = p1 * pow(2, 0) + p2 * pow(2, 1) + p3 * pow(2, 2);


    cout << "Nr blednego bitu "<< endl;
    cout << n << endl;

    if(n!=0)
        NEG(h, n-1);

    cout << endl << "Poprawiony strumien " << endl;

    for (int i = 0; i < 8; i++)
        cout << h[i] << endl;

    p4 = Parry(h, 7);    

    if (p4 != h[7])
    {
        cout << "Wykryto 2 bledy" << endl;
        return NULL;
    }

    return h;    
}


int main()
{   
    string data = "01000101";  
    
    string dlstr = data.substr(0, data.length() / 2);
    string drstr = data.substr(data.length() / 2);

    cout << "Dane pelnego bitu" << endl;
    cout << data.length() << " : " << data << endl;
    cout << "Lewy bit" << endl;
    cout << dlstr.length() << " : " << dlstr << endl;
    cout << "Prawy bit" << endl;
    cout << drstr.length() << " : " << drstr << endl << endl;


    cout << "Kod hamminga dla lewej czesci" << endl;
    int * k = CodeHamming(dlstr);


    cout << "Kod hamminga z bledem" << endl;

    NEG(k, 0);

    //NEG(k, 0);

    for (int i = 0; i < 7; i++)
        cout << k[i] << endl;
    cout << endl;

    DecodeHamming(k);

    cout << "=============================================" << endl;
    cout << "Kodowanie SECDED" << endl;


    k = CodeSECDED(dlstr);
    cout << endl;

    NEG(k, 1);
    NEG(k, 3);    

    cout << "Kod SECDED z bledem" << endl;
    for (int i = 0; i < 8; i++)
        cout << k[i] << endl;    
    cout << endl;

    DecodeSECDED(k);
}