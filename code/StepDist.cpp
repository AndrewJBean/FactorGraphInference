#include "Solvers.h"
#include <iostream>
#include <ctime>
#include <random>
using namespace std;

int const NumThreads = 16;
int const AlgoIters = 30;

namespace
{
    void MyPut(FloatType value)
    {
        fwrite(&(value), sizeof(value), 1, stdout);
    }
    void Normalize(vector<FloatType> & vec)
    {
        FloatType Total = 0.0;
        for(auto i:vec)
            Total+=i;
        for(auto & i:vec)
            i/=Total;
    }
};

int main(int argc, char** argv)
{
    cerr << "Read the factor graph from stdin." << endl;
    FactorGraph JustTheGraph;
    JustTheGraph.ReadGraphBin();

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    SumProductBP OneSolver;

    cerr << "Convergence of " << OneSolver.IDString() << "." << endl;
    cerr << OneSolver.PrintParams() << endl;

    vector< FloatType > theDist;
    theDist.resize(AlgoIters);

    clock_t begin_time = clock();
    long WallTime = time(NULL);

    // This includes any precomputation that is necessary to set up the algorithm
    //  this compute time will be included in the timestamps
    OneSolver.UseGraph(JustTheGraph);
    for(int i=0;i<AlgoIters;i++)
    {
        OneSolver.IterThreads( NumThreads );

        cerr << "                      \r"
            << i << "  "
            << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "   "
            << time(NULL)-WallTime ;

        OneSolver.GetDiffDist(theDist);
        for(auto val:theDist)
        {
            MyPut(val);
        }
        Normalize(theDist);
        for(auto val:theDist)
        {
            MyPut(val);
        }
    }
    cerr << endl;

    /////////////////////////////////////////////////////////////////////////////////

    return 0;
}

/*

name = 'ConesStepDist';
fid = fopen(name);
data = fread(fid,[1 inf],'*double')';
fclose(fid);
data = data(1:50*floor(numel(data)/50));
data = reshape(data,50,[]);
data1 = data(:,1:2:end);
data2 = data(:,2:2:end);
figure;semilogy(data1(:),'r');hold on;semilogy(data1(:),'.');
figure;semilogy(data2(:),'r');hold on;semilogy(data2(:),'.');

*/








