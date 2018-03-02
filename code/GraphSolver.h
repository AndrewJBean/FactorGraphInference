// This abstract base class defines the interface
//   that a graph solving algorithm should implement.


#pragma once
#include "FloatType.h"
#include <vector>
#include <string>
#include "DiscreteFactorGraph.h"
using namespace std;

class GraphSolver 
{
public:
	GraphSolver(){};																// DONE
	GraphSolver(FactorGraph newGraph);												// DONE
	virtual ~GraphSolver(){};														// DONE
	virtual string IDString()=0;													// DONE
	virtual string PrintParams()=0;													// DONE

	int  ReadGraphBin();															// DONE
	int  WriteGraphBin();															// DONE
	int  ReadGraphBin(string const &);												// DONE
	int  WriteGraphBin(string const &);												// DONE
	int  ReadGraph();																// DONE
	int  WriteGraph();																// DONE
	virtual int  UseGraph(FactorGraph newGraph);									// DONE
	virtual void ResetState()=0;													// DONE
	virtual void GetState(vector< vector< FloatType > > &)=0;						// DONE
	// RestoreState must first have the correct graph loaded
	virtual void RestoreState(vector< vector< FloatType > > const &)=0;				// DONE
	int  GetNumVars(){return graph.variables.size();}								// DONE
	int  GetNumFns(){return graph.functions.size();}								// DONE
	virtual void Iter()=0;															// DONE
	virtual void IterThreads(int NumThreads)=0;										// DONE

	virtual vector<FloatType> MarginalOfVar(int WhichVar)=0;						// DONE
	virtual void MarginalOfVars(vector< vector<FloatType> > & Result)=0;			// DONE
	virtual void DecisionOfVars(vector< int > & Result)=0;							// DONE

	int FullToMarginalIndex(int Fn, int Leave, int FullIndex);						// DONE
	int CoordsToMarginalIndex(int Fn, int Leave, vector< int > const & theCoords);	// DONE
	int CoordsToFullIndex(int Fn, vector< int > const & theCoords);					// DONE
	void IncSymbols(int Fn, vector< int > & theCoords);								// DONE
	void IncSymbolsExcept(int Fn, vector< int > & theCoords, int Except);			// DONE

	FactorGraph graph;
	vector< int > VarOfEdge;
	vector< int > FnOfEdge;

	vector< vector< int > > Var_Fn_Edge;
	vector< vector< int > > Fn_Var_Edge;
};
