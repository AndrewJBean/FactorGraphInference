// Select Beta to determine what fraction of the full update to make
// values closer to 1 for Beta indicate increased complexity, closer to Sum-Product
// subset size is adaptive
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

class BetaStepBP : public GraphSolver
{
public:
    BetaStepBP():Beta(0.5){};                                                       // DONE
    BetaStepBP(FactorGraph newGraph);                                               // DONE
    BetaStepBP(FloatType BetaParam);                                                // DONE
    BetaStepBP(FactorGraph newGraph,FloatType BetaParam);                           // DONE
    ~BetaStepBP(){};                                                                // DONE
    string IDString(){return string("BetaStepBP");};
    string PrintParams();

    int  UseGraph(FactorGraph newGraph);                                            // DONE
    void ResetState();                                                              // DONE
    void GetState(vector< vector< FloatType > > &);                                 // DONE
    void RestoreState(vector< vector< FloatType > > const &);                       // DONE
    void Iter();                                                                    // DONE
    struct ThreadData;                                                              // DONE
    static void * EdgeToFnThread(void *);                                           // DONE
    static void * EdgeToVarThread(void *);                                          // DONE
    static void * BeliefFromVarThread(void *);                                      // DONE
    void IterThreads(int NumThreads);                                               // DONE

    void UpdateMsgFromVar_Edge(int);                                                // DONE
    void UpdateMsgFromFn_Edge(int);                                                 // DONE
    void UpdateBeliefFromVar_Edge(int);                                             // DONE

    vector<FloatType> MarginalOfVar(int WhichVar);                                  // DONE
    void MarginalOfVars(vector< vector<FloatType> > & Result);                      // DONE
    void DecisionOfVars(vector< int > & Result);                                    // DONE

    vector<int> KUpdate;
    float Beta;

    // card is short for cardinality, i.e. size of the alphabet for that variable

    vector< vector< FloatType > >   EdgeToFn_Dist;          // NumEdges X (src var card)
    vector< vector< FloatType > >   EdgeToFn_DistEst;       // NumEdges X (src var card)
    vector< vector< FloatType > >   EdgeToFn_DistEst_old;   // NumEdges X (src var card)
    vector< vector< int > >         EdgeToFn_FindK;         // NumEdges X (src var card)

    vector< vector< FloatType > >   EdgeToVar_Marginal;     // NumEdges X (dest var card)
    vector< vector< FloatType > >   EdgeToVar_Dist;         // NumEdges X (dest var card)

    // temporary var for indexing into arrays
    vector< vector< int > >         EdgeToVar_Indices;      // NumEdges X (fn degree)

    struct myCompare
    {
        vector< FloatType > * vec;
        bool operator() (int i,int j) { return ((*vec)[i] > (*vec)[j]);}
    };
};


