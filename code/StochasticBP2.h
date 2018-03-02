// Modified (Simplified) Stochastic Belief Propagation, similar to Noorshams, Wainwright paper
// sampling is only applied on the variable to function messages
// This does not use the scalings as defined in the paper
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

class StochasticBP2 : public GraphSolver
{
public:
	StochasticBP2(){};																// DONE
	StochasticBP2(default_random_engine * RNG):r(RNG){};							// DONE
	StochasticBP2(FactorGraph newGraph);											// DONE
	~StochasticBP2(){};																// DONE
	string IDString(){return string("StochasticBP2");};
	string PrintParams();

	int  UseGraph(FactorGraph newGraph);											// DONE
	void ResetState();																// DONE
	void GetState(vector< vector< FloatType > > &);									// DONE
	void RestoreState(vector< vector< FloatType > > const &);						// DONE
	void Iter();																	// DONE
	void Iter(FloatType StepSize);													// DONE
	void IterThreads(int NumThreads);												// DONE
	void IterThreads(FloatType StepSize, int NumThreads);							// DONE

	void UpdateMsgFromVar_Edge(int);												// DONE
	void UpdateMsgFromFn_Edge(int,FloatType,
				default_random_engine *,uniform_real_distribution<FloatType> *);	// DONE

	vector<FloatType> MarginalOfVar(int WhichVar);									// DONE
	void MarginalOfVars(vector< vector<FloatType> > & Result);						// DONE
	void DecisionOfVars(vector< int > & Result);									// DONE

	// card is short for cardinality, i.e. size of the alphabet for that variable
	vector< vector< FloatType > > EdgeToFnMsg; // NumEdges X (src var card), same as the following
	vector< vector< FloatType > > EdgeToVarMsg; // NumEdges X (dest var card), same as the previous
	vector< vector< int > > EdgeToVarSymbols; // NumEdges X (fn degree)

	void SetRNG(default_random_engine * theRNG){r=theRNG;};							// DONE
	default_random_engine * r;
	int IterNumber;
};


