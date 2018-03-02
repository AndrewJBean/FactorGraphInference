// This is the standard Sum-Product Belief Propagation algorithm

#pragma once

#include "GraphSolver.h"
#include <vector>
#include <string>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class SumProductBP : public GraphSolver
{
public:
    SumProductBP(){};                                                               // DONE
    SumProductBP(FactorGraph newGraph){UseGraph(newGraph);};                        // DONE
    ~SumProductBP(){};                                                              // DONE
    string IDString(){return string("SumProductBP");};
    string PrintParams();

    int  UseGraph(FactorGraph newGraph);                                            // DONE
    void ResetState();                                                              // DONE
    void GetState(vector< vector< FloatType > > &);                                 // DONE
    void RestoreState(vector< vector< FloatType > > const &);                       // DONE
    void Iter();                                                                    // DONE
    struct ThreadData;                                                              // DONE
    static void * EdgeToFnThread(void * d);                                         // DONE
    static void * EdgeToVarThread(void * d);                                        // DONE
    void IterThreads(int NumThreads);                                               // DONE

    void UpdateMsgFromVar_Edge(int);                                                // DONE
    void UpdateMsgFromFn_Edge(int);                                                 // DONE

    vector<FloatType> MarginalOfVar(int WhichVar);                                  // DONE
    void MarginalOfVars(vector< vector<FloatType> > & Result);                      // DONE
    void DecisionOfVars(vector< int > & Result);                                    // DONE

    void GetDiffDist(vector<FloatType> &);                                          //

    // card is short for cardinality, i.e. size of the alphabet for that variable
    vector< vector< FloatType > > EdgeToFnMsg;         // NumEdges X (src var card), same as the following
    vector< vector< FloatType > > EdgeToFnMsg_old;     // NumEdges X (src var card), same as the following
    vector< vector< FloatType > > EdgeToVarMsg;        // NumEdges X (dest var card), same as the previous
    vector< vector< int > >       EdgeToVarSymbols;    // NumEdges X (fn degree)
    vector<FloatType> TempVec;
};


