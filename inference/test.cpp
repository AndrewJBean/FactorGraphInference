#include <iostream>
using namespace std;

int const w=240;
int const h=120;
int const maxval=30-1;

void ASCII()
{
    cout << "P2" << endl;
    cout << w << " " << h << endl ;
    cout << maxval << endl ;
    for(int i=0;i<h;i++)
    {
        for(int j=0;j<w;j++)
        {
            cout << (i-j + (maxval+1)*100)%(maxval+1) << " " ;
        }
        cout << endl;
    }
}

void BINARY()
{
    cout << "P5" << endl;
    cout << w << " " << h << endl ;
    cout << maxval << endl ;
    for(int i=0;i<h;i++)
    {
        for(int j=0;j<w;j++)
        {
            cout.put( char( (i-j+100*(maxval+1)) % (maxval+1) ) );
        }
    }
}

int main(int argc, char** argv)
{
    BINARY();
    return 0;
}

