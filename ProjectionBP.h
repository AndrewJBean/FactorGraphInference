// Select KUpdate, the number of elements to update for the var-to-fn msg estimates
// values are 1 to variable cardinality
// subset size is not adaptive
// compression is only applied on the variable to function messages
// the estimates of the variable to function messages are not enforced to sum to 1

#pragma once

#include "GraphSolver.h"
#include <vector>
#include <string>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class ProjectionBP : public GraphSolver
{
public:
	ProjectionBP():KUpdate(1){};														// DONE
	ProjectionBP(FactorGraph newGraph):KUpdate(1){UseGraph(newGraph);};					// DONE
	ProjectionBP(int KParam):KUpdate(KParam){};											// DONE
	ProjectionBP(FactorGraph newGraph,int KParam):KUpdate(KParam){UseGraph(newGraph);};	// DONE
	~ProjectionBP(){};																	// DONE
	string IDString(){return string("ProjectionBP");};
	string PrintParams();

	int  UseGraph(FactorGraph newGraph);												// DONE
	void ResetState();																	// DONE
	void GetState(vector< vector< FloatType > > &);										// 
	void RestoreState(vector< vector< FloatType > > const &);							// 
	void Iter();																		// DONE
	struct ThreadData;																	// DONE
	static void * EdgeToFnThread(void *);												// DONE
	static void * EdgeToVarThread(void *);												// DONE
	static void * BeliefFromVarThread(void *);											// DONE
	void IterThreads(int NumThreads);													// DONE

	void UpdateMsgFromVar_Edge(int);													// DONE
	void UpdateMsgFromFn_Edge(int);														// DONE
	void UpdateBeliefFromVar_Edge(int);													// DONE

	vector<FloatType> MarginalOfVar(int WhichVar);										// DONE
	void MarginalOfVars(vector< vector<FloatType> > & Result);							// DONE
	void DecisionOfVars(vector< int > & Result);										// DONE

	int partition_largest(		vector< FloatType > const & input,
								vector< int > & Scratch,
								int p, int r);											// DONE
	void quick_select_largest(	vector< FloatType > const & input,
								vector< int > & Scratch,
								int p, int r, int k);									// DONE

	int KUpdate;
	
	// card is short for cardinality, i.e. size of the alphabet for that variable

	vector< vector< FloatType > >	EdgeToFn_Dist;				// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_DistEst;			// NumEdges X (src var card)
	vector< vector< FloatType > >	EdgeToFn_DistEst_old;		// NumEdges X (src var card)
	vector< vector< int > >			EdgeToFn_FindK;				// NumEdges X (src var card)

	vector< vector< FloatType > >	EdgeToVar_Marginal;			// NumEdges X (dest var card)
	vector< vector< FloatType > >	EdgeToVar_Dist;				// NumEdges X (dest var card)

	// temporary var for indexing into arrays
	vector< vector< int > > 		EdgeToVar_Indices;			// NumEdges X (fn degree)
};


