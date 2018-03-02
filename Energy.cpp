#include "Solvers.h"
#include <iostream>
#include <ctime>
#include <random>
using namespace std;

int const NumThreads = 1;
int const AlgoIters = 1000;

int main(int argc, char** argv)
{
	cerr << "Read the factor graph from stdin." << endl;
	FactorGraph JustTheGraph;
	JustTheGraph.ReadGraphBin();

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
	mySolvers.push_back(new SumProductBP());				// 1
	// for(auto p:BetaParam)									// 2:17
	// 	mySolvers.push_back(new BetaStepBP(p));
	for(auto p:KParam)										// 18:27
		mySolvers.push_back(new ProjectionBP(p));
	// for(auto p:KParam)										// 28:37
	// 	mySolvers.push_back(new AProjectionBP(p));
	// for(auto p:KParam)										// 38:47
	// 	mySolvers.push_back(new KProjectionBP(p));
	// mySolvers.push_back(new ZoomBP());						// 48
	mySolvers.push_back(new FastZoomBP());					// 49
	mySolvers.push_back(new StochasticBP(&generator));		// 50
	mySolvers.push_back(new StochasticBP2(&generator));		// 51


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
		for(int i=0;i<AlgoIters;i++)
		{
			OneSolver->IterThreads( NumThreads );

			cerr << "                      \r"
				<< i << "  " 
				<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "   " 
				<< time(NULL)-WallTime ;

			OneSolver->DecisionOfVars(GraphDecisions);

			cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "," ;
			cout << time(NULL)-WallTime << "," ;
			cout << -JustTheGraph.LogValue(GraphDecisions) << endl;
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










