#include "Solvers.h"
#include <iostream>
#include <ctime>
#include <random>
using namespace std;

int const NumThreads = 16;
int const AlgoIters = 1000;
float MaxRuntime = 18000; // seconds

//
// Usage:
//  ./Disparities inputfile
//  OR
//  ./Disparities inputfile /output/path/
//  OR
//  ./Disparities < inputfile
//
// inputfile should contain data for the factor graph plus 3 more integers: w, h, and D
// output images go to ./images/ folder or /output/path/images/ folder
//

int main(int argc, char** argv)
{
    cerr << "Read the factor graph." << endl;
    FactorGraph JustTheGraph;
    int width,height,D;

    if(argc>1)
    {
        string FileName(argv[1]);
        ifstream file (FileName.c_str(), std::ifstream::binary);
        file.seekg(0);
        JustTheGraph.ReadGraphBin(file);
        FactorGraph::MyGet(width,file);
        FactorGraph::MyGet(height,file);
        FactorGraph::MyGet(D,file);
        file.close();
    }
    else
    {
        JustTheGraph.ReadGraphBin();
        FactorGraph::MyGet(width);
        FactorGraph::MyGet(height);
        FactorGraph::MyGet(D);
    }

    // vector<FloatType> BetaParam = {0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95};
    // vector<int> KParam = {1,2,3,4,5,6,7,8,9,10};
    vector<FloatType> BetaParam = {0.95};
    vector<int> KParam = {1};
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    default_random_engine generator;
    vector<GraphSolver*> mySolvers ;
    mySolvers.push_back(new SumProductBP());
    // for(auto p:BetaParam)
    //  mySolvers.push_back(new BetaStepBP(p));
    for(auto p:KParam)
        mySolvers.push_back(new ProjectionBP(p));
    // for(auto p:KParam)
    //  mySolvers.push_back(new AProjectionBP(p));
    // for(auto p:KParam)
    //  mySolvers.push_back(new KProjectionBP(p));
    // mySolvers.push_back(new ZoomBP());
    mySolvers.push_back(new FastZoomBP());
    mySolvers.push_back(new StochasticBP(&generator));
    mySolvers.push_back(new StochasticBP2(&generator));

    string RootPath = "";
    if(argc>2)
        RootPath = argv[2];
    cerr << "RootPath=" << RootPath << endl;


    for(auto OneSolver:mySolvers)
    {
        cerr << "Convergence of " << OneSolver->IDString() << "." << endl;
        cerr << OneSolver->PrintParams() << endl;

        vector<int> GraphDecisions;
        clock_t begin_time = clock();
        long WallTime = time(NULL);

        // This includes any precomputation that is necessary to set up the algorithm
        //  this compute time will be included in the timestamps
        OneSolver->UseGraph(JustTheGraph);
        for(int i=0;i<AlgoIters && (float( clock () - begin_time ) /  CLOCKS_PER_SEC)<MaxRuntime ;i++)
        {
            OneSolver->IterThreads( NumThreads );

            cerr << "\r"
                << i << "  "
                << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "   "
                << time(NULL)-WallTime << "    " ;

            OneSolver->DecisionOfVars(GraphDecisions);

            string Number = to_string(i+1);
            if(i+1 < 1000) Number = string("0") + Number;
            if(i+1 < 100) Number = string("0") + Number;
            if(i+1 < 10) Number = string("0") + Number;
            string FileName = RootPath+"images/" + OneSolver->IDString() + "-" + OneSolver->PrintParams() + "_" + Number + ".pgm" ;
            // string FileName = OneSolver->IDString() + "-" + OneSolver->PrintParams() + ".pgm" ;
            ofstream file (FileName.c_str(), ios::out|ios::binary|ios::trunc);
            if(file.is_open())
            {
                file.seekp (0);
                file << "P5" << endl;
                file << width << " " << height << endl ;
                file << D-1 << endl ;
                for(int looph=0;looph<height;looph++)
                {
                    for(int loopw=0;loopw<width;loopw++)
                    {
                        file.put( char( D-1-GraphDecisions[loopw + width*(height-looph-1)] ) );
                    }
                }
                file.close();
            }
            else
                cerr << "SOMETHING WRONG: " << FileName << endl ;

            // cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "," ;
            // cout << time(NULL)-WallTime << "," ;
            // cout << -JustTheGraph.LogValue(GraphDecisions) << endl;
        }
        cerr << endl;
        delete OneSolver;
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    return 0;
}










