// Only one element is updated for the msg estimates
// compression is applied to var-to-fn messages and fn-to-var messages
// the estimates of the variable to function messages ARE NOT enforced to sum to 1
// Simply choose the message element that is the most off from the estimate

#pragma once

#include "GraphSolver.h"
#include <vector>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class NewZoomBP : public GraphSolver
{
public:
	NewZoomBP():QuantFnToVar(1),QuantVarToFn(1){};											// DONE
	NewZoomBP(FactorGraph newGraph):QuantFnToVar(1),QuantVarToFn(1){UseGraph(newGraph);};	// DONE
	~NewZoomBP(){};																			// DONE
	string IDString(){return string("NewZoomBP");};
	string PrintParams();

	int  UseGraph(FactorGraph newGraph);													// DONE
	void ResetState();																		// DONE
	void GetState(vector< vector< FloatType > > &);											// DONE
	void RestoreState(vector< vector< FloatType > > const &);								// DONE
	void Iter();																			// DONE
	struct ThreadData;																		// DONE
	static void * EdgeToFnThread(void *);													// DONE
	static void * EdgeToVarThread(void *);													// DONE
	static void * BeliefFromVarThread(void *);												// DONE
	void IterThreads(int NumThreads);														// DONE

	void UpdateMsgFromVar_Edge(int);														// DONE
	void UpdateMsgFromFn_Edge(int);															// DONE
	void UpdateBeliefFromVar_Edge(int);														// DONE

	vector<FloatType> MarginalOfVar(int WhichVar);											// DONE
	void MarginalOfVars(vector< vector<FloatType> > & Result);								// DONE
	void DecisionOfVars(vector< int > & Result);											// DONE

	// card is short for cardinality, i.e. size of the alphabet for that variable

	vector< vector< FloatType > >	EdgeToFn_NonNormal;		// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_Dist;			// NumEdges X (src var card)
	vector< int >					EdgeToFn_I;
	vector< FloatType >				EdgeToFn_D;
	vector< int >					EdgeToFn_Q;
	vector< FloatType >				EdgeToFn_R;
	vector< FloatType >				EdgeToFn_Alpha;
	vector< FloatType >				EdgeToFn_Gamma;
	vector< vector< FloatType > >	EdgeToFn_DistEst;		// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_DistEst_old;	// NumEdges X (src var card)

	vector< vector< FloatType > >	EdgeToVar_Marginal;		// NumEdges X (dest var card)
	vector< vector< FloatType > >	EdgeToVar_Dist;			// NumEdges X (dest var card)
	vector< int >					EdgeToVar_I;
	vector< FloatType >				EdgeToVar_D;
	vector< int >					EdgeToVar_Q;
	vector< FloatType >				EdgeToVar_R;
	vector< FloatType >				EdgeToVar_Alpha;
	vector< FloatType >				EdgeToVar_Gamma;
	vector< vector< FloatType > >	EdgeToVar_DistEst;		// NumEdges X (dest var card)
	vector< vector< FloatType > >	EdgeToVar_DistEst_old;	// NumEdges X (dest var card)

	// temporary var for indexing into arrays
	vector< vector< int > > 		EdgeToVar_Indices;		// NumEdges X (fn degree)
};


