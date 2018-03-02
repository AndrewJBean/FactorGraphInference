// Stochastic Belief Propagation, as in the Noorshams, Wainwright paper
// sampling is only applied on the variable to function messages
// This uses the scalings as defined in the paper
// their algorithm has been generalized for higher order fn nodes

#pragma once

#include "GraphSolver.h"
#include <vector>
#include <random>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class StochasticBP : public GraphSolver
{
public:
	StochasticBP(){};																// DONE
	StochasticBP(default_random_engine * RNG):r(RNG){};								// DONE
	StochasticBP(FactorGraph newGraph);												// DONE
	~StochasticBP(){};																// DONE
	string IDString(){return string("StochasticBP");};
	string PrintParams();

	int  UseGraph(FactorGraph newGraph);											// DONE
	void ResetState();																// DONE
	void GetState(vector< vector< FloatType > > &);									// DONE
	void RestoreState(vector< vector< FloatType > > const &);						// DONE
	void Precompute(int Edge);														// DONE
	void Iter();																	// DONE
	void Iter(FloatType StepSize);													// DONE
	struct ThreadData;																// DONE
	static void * EdgeToFnThread(void * d);											// DONE
	static void * EdgeToVarThread(void * d);										// DONE
	void IterThreads(int NumThreads);												// DONE
	void IterThreads(FloatType StepSize, int NumThreads);							// DONE

	void UpdateMsgFromVar_Edge(int);												// DONE
	void UpdateMsgFromFn_Edge(int,FloatType,
				default_random_engine *,uniform_real_distribution<FloatType> *);	// DONE

	vector<FloatType> MarginalOfVar(int WhichVar);									// DONE
	void MarginalOfVars(vector< vector<FloatType> > & Result);						// DONE
	void DecisionOfVars(vector< int > & Result);									// DONE

	// card is short for cardinality, i.e. size of the alphabet for that variable
	vector< vector< FloatType > > Betas; // NumEdges X (Prod card / var card)
	vector< FloatType > BetaMax; // NumEdges X 1
	vector< vector< FloatType > > Gammas; // NumEdges X (Prod card)
	vector< vector< FloatType > > EdgeToFnMsg; // NumEdges X (src var card), same as the following
	vector< vector< FloatType > > EdgeToVarMsg; // NumEdges X (dest var card), same as the previous
	vector< vector< int > > EdgeToVarSymbols; // NumEdges X (fn degree)

	void SetRNG(default_random_engine * theRNG){r=theRNG;};							// DONE

	int IterNumber;
	default_random_engine * r;
};


