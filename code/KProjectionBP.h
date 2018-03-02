// Select KUpdate, the number of elements to update for the msg estimates
// values are 1 to variable cardinality
// subset size is not adaptive
// compression is applied to var-to-fn and fn-to-var messages
// the estimates of the variable to function messages ARE enforced to sum to 1

#pragma once

#include "GraphSolver.h"
#include <vector>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class KProjectionBP : public GraphSolver
{
public:
	KProjectionBP():KUpdate(1){};															// DONE
	KProjectionBP(FactorGraph newGraph):KUpdate(1){UseGraph(newGraph);};					// DONE
	KProjectionBP(int KParam):KUpdate(KParam){};											// DONE
	KProjectionBP(FactorGraph newGraph,int KParam):KUpdate(KParam){UseGraph(newGraph);};	// DONE									// DONE
	~KProjectionBP(){};																		// DONE
	string IDString(){return string("KProjectionBP");};
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

	int partition_largest(		vector< FloatType > const & input,
								vector< int > & Scratch,
								int p, int r);												// DONE
	void quick_select_largest(	vector< FloatType > const & input,
								vector< int > & Scratch,
								int p, int r, int k);										// DONE

	int KUpdate;
	
	// card is short for cardinality, i.e. size of the alphabet for that variable

	vector< vector< FloatType > >	EdgeToFn_NonNormal;					// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_Dist;						// NumEdges X (src var card)
	vector< FloatType >				EdgeToFn_Alpha;
	vector< FloatType >				EdgeToFn_Gamma;
	vector< vector< FloatType > >	EdgeToFn_DistEst;					// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_DistEst_old;				// NumEdges X (src var card)
	vector< vector< int > >			EdgeToFn_FindK;						// NumEdges X (src var card)

	vector< vector< FloatType > >	EdgeToVar_Marginal;					// NumEdges X (dest var card)
	vector< vector< FloatType > >	EdgeToVar_Dist;						// NumEdges X (dest var card)
	vector< FloatType >				EdgeToVar_Alpha;
	vector< FloatType >				EdgeToVar_Gamma;
	vector< vector< FloatType > >	EdgeToVar_DistEst;					// NumEdges X (dest var card)
	vector< vector< FloatType > >	EdgeToVar_DistEst_old;				// NumEdges X (dest var card)
	vector< vector< int > >			EdgeToVar_FindK;					// NumEdges X (dest var card)

	// temporary var for indexing into arrays
	vector< vector< int > > 		EdgeToVar_Indices;					// NumEdges X (fn degree)
};


