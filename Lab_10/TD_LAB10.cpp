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
        {
            matrix[i] = new int[y];
        }
    }

    Matrix(int x, int y, string data)
    {
        this->x = x;
        this->y = y;

        if (x * y != data.length())        
            throw ("Data length error");            
        
        matrix = new int* [x];
        for (int i = 0; i < x; i++)
        {
            matrix[i] = new int [y];
        }

        int count = 0;
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
            {
                matrix[i][j] = data[count++] - '0';
            }            
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
        {
            matrix[i] = new int[y];
        }

        int count = 0;
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
            {
                matrix[i][j] = data[count++];
            }
        }
    }

    void show()
    {
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
            {
                cout << matrix[i][j]<< " ";
            }
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
                    {
                        suma = suma + this->matrix[i][h] * b.matrix[h][j];
                    }
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
};


void NEG(int* data,int number)
{
    data[number] = !data[number];
}

int main()
{
    
    string gstr = "1101101110000111010000100001";
    string hstr = "101010101100110001111";

    string dstr = "1101"; 
 

    Matrix G(7, 4, gstr);    
    Matrix H(3, 7, hstr);

    Matrix D(4, 1, dstr);
     
    D.show(); 
        
    Matrix K = G * D;    
    K % 2;
    K.show();

    int* kp = K.list();
    NEG(kp,1);

    Matrix Kp(7, 1, kp, 7);
    Kp.show();

    Matrix S = H * Kp;
    S % 2;
    S.show();

}