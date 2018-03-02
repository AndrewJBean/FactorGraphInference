#include "SumProductBP.h"
#include <cmath>
#include <pthread.h>
#include <algorithm>

#ifndef NULL
#define NULL 0
#endif

int const ChunkSize = 50;

string SumProductBP::PrintParams()
{
    return string("none");
}

//  __    __                       ______                                 __          ___  ___
// |  \  |  \                     /      \                               |  \        /   \|   \
// | $$  | $$  _______   ______  |  $$$$$$\  ______    ______    ______  | $$____   /  $$$ \$$$\
// | $$  | $$ /       \ /      \ | $$ __\$$ /      \  |      \  /      \ | $$    \ |  $$     \$$\
// | $$  | $$|  $$$$$$$|  $$$$$$\| $$|    \|  $$$$$$\  \$$$$$$\|  $$$$$$\| $$$$$$$\| $$      | $$
// | $$  | $$ \$$    \ | $$    $$| $$ \$$$$| $$   \$$ /      $$| $$  | $$| $$  | $$| $$      | $$
// | $$__/ $$ _\$$$$$$\| $$$$$$$$| $$__| $$| $$      |  $$$$$$$| $$__/ $$| $$  | $$ \$$\_  _/  $$
//  \$$    $$|       $$ \$$     \ \$$    $$| $$       \$$    $$| $$    $$| $$  | $$  \$$ \|   $$
//   \$$$$$$  \$$$$$$$   \$$$$$$$  \$$$$$$  \$$        \$$$$$$$| $$$$$$$  \$$   \$$   \$$$ \$$$
//                                                             | $$
//                                                             | $$
//                                                              \$$



// DONE
// Verified by inspection only.
int SumProductBP::UseGraph(FactorGraph newGraph)
{
    GraphSolver::UseGraph(newGraph);

    int EdgeCount = VarOfEdge.size();

    EdgeToVarMsg.resize(EdgeCount);
    EdgeToFnMsg.resize(EdgeCount);
    EdgeToVarSymbols.resize(EdgeCount);

    for(int WhichEdge=0;WhichEdge<EdgeCount;WhichEdge++)
    {
        int tempSize = graph.variables[ VarOfEdge[WhichEdge] ].cardinality;
        EdgeToFnMsg[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) );
        EdgeToVarMsg[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) );
        EdgeToVarSymbols[WhichEdge].resize(graph.functions[FnOfEdge[WhichEdge]].vars.size());
    }
    EdgeToFnMsg_old = EdgeToFnMsg;
    return 0;
}





void SumProductBP::ResetState()
{
    for(int WhichEdge=0;WhichEdge<EdgeToFnMsg.size();WhichEdge++)
    {
        int tempSize = graph.variables[ VarOfEdge[WhichEdge] ].cardinality;
        EdgeToVarMsg[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) );
    }
    EdgeToFnMsg_old = EdgeToVarMsg;
}



void SumProductBP::GetState(vector< vector< FloatType > > & State)
{
    // State consists of:
    //  - EdgeToVarMsg[][]
    State.resize(EdgeToVarMsg.size());
    int Index=0;
    for(auto & v:EdgeToVarMsg)
    {
        State[Index++] = v;
    }
}


void SumProductBP::RestoreState(vector< vector< FloatType > > const & State)
{
    int Index=0;
    for(auto & v:EdgeToVarMsg)
    {
        v = State[Index++];
    }
}





//  ______    __                                      __      __
// |      \  |  \                                    |  \    |  \
//  \$$$$$$ _| $$_     ______    ______    ______   _| $$_    \$$  ______   _______
//   | $$  |   $$ \   /      \  /      \  |      \ |   $$ \  |  \ /      \ |       \
//   | $$   \$$$$$$  |  $$$$$$\|  $$$$$$\  \$$$$$$\ \$$$$$$  | $$|  $$$$$$\| $$$$$$$\
//   | $$    | $$ __ | $$    $$| $$   \$$ /      $$  | $$ __ | $$| $$  | $$| $$  | $$
//  _| $$_   | $$|  \| $$$$$$$$| $$      |  $$$$$$$  | $$|  \| $$| $$__/ $$| $$  | $$
// |   $$ \   \$$  $$ \$$     \| $$       \$$    $$   \$$  $$| $$ \$$    $$| $$  | $$
//  \$$$$$$    \$$$$   \$$$$$$$ \$$        \$$$$$$$    \$$$$  \$$  \$$$$$$  \$$   \$$



//  ________                                  __      __
// |        \                                |  \    |  \
// | $$$$$$$$ __    __  _______    _______  _| $$_    \$$  ______   _______    _______
// | $$__    |  \  |  \|       \  /       \|   $$ \  |  \ /      \ |       \  /       \
// | $$  \   | $$  | $$| $$$$$$$\|  $$$$$$$ \$$$$$$  | $$|  $$$$$$\| $$$$$$$\|  $$$$$$$
// | $$$$$   | $$  | $$| $$  | $$| $$        | $$ __ | $$| $$  | $$| $$  | $$ \$$    \
// | $$      | $$__/ $$| $$  | $$| $$_____   | $$|  \| $$| $$__/ $$| $$  | $$ _\$$$$$$\
// | $$       \$$    $$| $$  | $$ \$$     \   \$$  $$| $$ \$$    $$| $$  | $$|       $$
//  \$$        \$$$$$$  \$$   \$$  \$$$$$$$    \$$$$  \$$  \$$$$$$  \$$   \$$ \$$$$$$$




void SumProductBP::Iter()
{
    int const NumMessages = EdgeToFnMsg.size();

    // update every EdgeToFnMsg
    for(int i=0;i<NumMessages;i++)
    {
        UpdateMsgFromVar_Edge(i);
    }

    // update every EdgeToVarMsg
    for(int i=0;i<NumMessages;i++)
    {
        UpdateMsgFromFn_Edge(i);
    }
}





///////////////////////
struct SumProductBP::ThreadData
{
    SumProductBP * Graph;
    int Progress;
    pthread_mutex_t Mutex;
};
///////////////////////

void * SumProductBP::EdgeToFnThread(void * d)
{
    // cast the input pointer so we can use it
    ThreadData * theData = (ThreadData *)d;

    int QuitThread;
    do
    {
        int StartCompute;
        QuitThread = 0;
        pthread_mutex_lock(&(theData->Mutex));
        if(theData->Progress < theData->Graph->EdgeToFnMsg.size())
        {
            StartCompute = theData->Progress;
            theData->Progress += ChunkSize;
        }
        else QuitThread = 1;
        pthread_mutex_unlock(&(theData->Mutex));
        if(!QuitThread)
        {
            for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToFnMsg.size();i++)
                theData->Graph->UpdateMsgFromVar_Edge(i);
        }
    }while(!QuitThread);
    return NULL;
}


void * SumProductBP::EdgeToVarThread(void * d)
{
    // cast the input pointer so we can use it
    ThreadData * theData = (ThreadData *)d;

    int QuitThread;
    do
    {
        int StartCompute;
        QuitThread = 0;
        pthread_mutex_lock(&(theData->Mutex));
        if(theData->Progress < theData->Graph->EdgeToVarMsg.size())
        {
            StartCompute = theData->Progress;
            theData->Progress += ChunkSize;
        }
        else QuitThread = 1;
        pthread_mutex_unlock(&(theData->Mutex));
        if(!QuitThread)
        {
            for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToVarMsg.size();i++)
                theData->Graph->UpdateMsgFromFn_Edge(i);
        }
    }while(!QuitThread);
    return NULL;
}




void SumProductBP::IterThreads(int NumThreads)
{
    if(NumThreads==1)
    {
        Iter();
        return;
    }
    //////////////////////////////////////////////
    // update every EdgeToFnMsg
    {
        vector<pthread_t> myThreads;
        myThreads.resize(NumThreads);

        ThreadData theData;

        theData.Graph = this;
        theData.Progress = 0;
        pthread_mutex_init(&theData.Mutex,NULL);

        for(int i=0;i<NumThreads;i++)
        {
            pthread_create( &myThreads[i], NULL, EdgeToFnThread, (void*)(&theData) );
        }
        for(int i=0;i<NumThreads;i++)
        {
            pthread_join(myThreads[i],NULL);
        }
        pthread_mutex_destroy(&theData.Mutex);
    }
    //////////////////////////////////////////////
    // update every EdgeToVarMsg
    {
        vector<pthread_t> myThreads;
        myThreads.resize(NumThreads);

        ThreadData theData;

        theData.Graph = this;
        theData.Progress = 0;
        pthread_mutex_init(&theData.Mutex,NULL);

        for(int i=0;i<NumThreads;i++)
        {
            pthread_create( &myThreads[i], NULL, EdgeToVarThread, (void*)(&theData) );
        }
        for(int i=0;i<NumThreads;i++)
        {
            pthread_join(myThreads[i],NULL);
        }
        pthread_mutex_destroy(&theData.Mutex);
    }
    //////////////////////////////////////////////
}









//  ________        __                            __    __                  __              __
// |        \      |  \                          |  \  |  \                |  \            |  \
// | $$$$$$$$  ____| $$  ______    ______        | $$  | $$  ______    ____| $$  ______   _| $$_     ______
// | $$__     /      $$ /      \  /      \       | $$  | $$ /      \  /      $$ |      \ |   $$ \   /      \
// | $$  \   |  $$$$$$$|  $$$$$$\|  $$$$$$\      | $$  | $$|  $$$$$$\|  $$$$$$$  \$$$$$$\ \$$$$$$  |  $$$$$$\
// | $$$$$   | $$  | $$| $$  | $$| $$    $$      | $$  | $$| $$  | $$| $$  | $$ /      $$  | $$ __ | $$    $$
// | $$_____ | $$__| $$| $$__| $$| $$$$$$$$      | $$__/ $$| $$__/ $$| $$__| $$|  $$$$$$$  | $$|  \| $$$$$$$$
// | $$     \ \$$    $$ \$$    $$ \$$     \       \$$    $$| $$    $$ \$$    $$ \$$    $$   \$$  $$ \$$     \
//  \$$$$$$$$  \$$$$$$$ _\$$$$$$$  \$$$$$$$        \$$$$$$ | $$$$$$$   \$$$$$$$  \$$$$$$$    \$$$$   \$$$$$$$
//                     |  \__| $$                          | $$
//                      \$$    $$                          | $$
//                       \$$$$$$                            \$$
//  ________                                  __      __
// |        \                                |  \    |  \
// | $$$$$$$$ __    __  _______    _______  _| $$_    \$$  ______   _______    _______
// | $$__    |  \  |  \|       \  /       \|   $$ \  |  \ /      \ |       \  /       \
// | $$  \   | $$  | $$| $$$$$$$\|  $$$$$$$ \$$$$$$  | $$|  $$$$$$\| $$$$$$$\|  $$$$$$$
// | $$$$$   | $$  | $$| $$  | $$| $$        | $$ __ | $$| $$  | $$| $$  | $$ \$$    \
// | $$      | $$__/ $$| $$  | $$| $$_____   | $$|  \| $$| $$__/ $$| $$  | $$ _\$$$$$$\
// | $$       \$$    $$| $$  | $$ \$$     \   \$$  $$| $$ \$$    $$| $$  | $$|       $$
//  \$$        \$$$$$$  \$$   \$$  \$$$$$$$    \$$$$  \$$  \$$$$$$  \$$   \$$ \$$$$$$$


void SumProductBP::UpdateMsgFromVar_Edge(int i)
{
    EdgeToFnMsg_old[i] = EdgeToFnMsg[i];
    // "clear" the message to all ones
    for(int k=0;k<graph.variables[ VarOfEdge[i] ].cardinality;k++)
        EdgeToFnMsg[ i ] [ k ] = 1;
    // find edge to var messages attached to VarOfEdge[i]
    // go through all the edges/functions connected to the var
    for(int j=0;j<Var_Fn_Edge[ VarOfEdge[i] ].size();j++)
    {
        // one of them should be the current edge, don't use in update
        if(Var_Fn_Edge[ VarOfEdge[i] ][j] != i)
            for(int k=0;k<graph.variables[ VarOfEdge[i] ].cardinality;k++)
                EdgeToFnMsg[ i ] [k] *=
                    EdgeToVarMsg[ Var_Fn_Edge[ VarOfEdge[i] ][j] ] [k];
    }
}



void SumProductBP::UpdateMsgFromFn_Edge(int i)
{
    // which function does the edge come from?
    int theFn=FnOfEdge[i];

    // Find which neighbor variable to leave out of summations
    int LeaveOut = -10000; // likely to cause segfault if there's some inconsistency/bug
    for(int j=0;j<graph.functions[theFn].vars.size();j++)
    {
        if(Fn_Var_Edge[ theFn ][j] == i)
            LeaveOut=j;
    }

    // zero out the message
    for(int j=0;j<EdgeToVarMsg[i].size();j++)
    {
        EdgeToVarMsg[i][j] = 0.0;
    }

    // zero out the symbols
    for(int j=0;j<EdgeToVarSymbols[i].size();j++)
    {
        EdgeToVarSymbols[i][j] = 0;
    }

    // go through every function value
    for(int j=0;j<graph.functions[theFn].ValuesSize;j++)
    {
        FloatType theProduct = 1;
        // product of the incoming messages at the current point (exclude current var)
        for(int k=0;k<EdgeToVarSymbols[i].size();k++)
        {
            // EdgeToFnMsg [desired edge] [value of the variable]
            // desired edge is k'th one attached to the current
            //     function, hence Fn_Var_Edge[theFn][k]
            // value of the variable is a coordinate 'k' from EdgeToVarSymbols[i]
            // note: the vars attached to the fn and the var symbols are same order
            if(k != LeaveOut) theProduct *=
                EdgeToFnMsg[  Fn_Var_Edge[theFn][k]  ][  EdgeToVarSymbols[i][k]  ];
        }
        EdgeToVarMsg[i][ EdgeToVarSymbols[i][LeaveOut] ]
            += graph.values[graph.functions[theFn].values][j] * theProduct;
        IncSymbols(theFn,EdgeToVarSymbols[i]);
    }

    // normalize the EdgeToVarMsg
    FloatType MsgTotal = 0.0;
    for(int j=0;j<EdgeToVarMsg[i].size();j++)
    {
        MsgTotal += EdgeToVarMsg[i][j];
    }
    MsgTotal = 1.0 / MsgTotal;
    for(int j=0;j<EdgeToVarMsg[i].size();j++)
    {
        EdgeToVarMsg[i][j] *= MsgTotal;
    }
}













//  __    __    __      __  __  __    __
// |  \  |  \  |  \    |  \|  \|  \  |  \
// | $$  | $$ _| $$_    \$$| $$ \$$ _| $$_    __    __
// | $$  | $$|   $$ \  |  \| $$|  \|   $$ \  |  \  |  \
// | $$  | $$ \$$$$$$  | $$| $$| $$ \$$$$$$  | $$  | $$
// | $$  | $$  | $$ __ | $$| $$| $$  | $$ __ | $$  | $$
// | $$__/ $$  | $$|  \| $$| $$| $$  | $$|  \| $$__/ $$
//  \$$    $$   \$$  $$| $$| $$| $$   \$$  $$ \$$    $$
//   \$$$$$$     \$$$$  \$$ \$$ \$$    \$$$$  _\$$$$$$$
//                                           |  \__| $$
//                                            \$$    $$
//                                             \$$$$$$
//  ________                                  __      __
// |        \                                |  \    |  \
// | $$$$$$$$ __    __  _______    _______  _| $$_    \$$  ______   _______    _______
// | $$__    |  \  |  \|       \  /       \|   $$ \  |  \ /      \ |       \  /       \
// | $$  \   | $$  | $$| $$$$$$$\|  $$$$$$$ \$$$$$$  | $$|  $$$$$$\| $$$$$$$\|  $$$$$$$
// | $$$$$   | $$  | $$| $$  | $$| $$        | $$ __ | $$| $$  | $$| $$  | $$ \$$    \
// | $$      | $$__/ $$| $$  | $$| $$_____   | $$|  \| $$| $$__/ $$| $$  | $$ _\$$$$$$\
// | $$       \$$    $$| $$  | $$ \$$     \   \$$  $$| $$ \$$    $$| $$  | $$|       $$
//  \$$        \$$$$$$  \$$   \$$  \$$$$$$$    \$$$$  \$$  \$$$$$$  \$$   \$$ \$$$$$$$




// DONE
vector<FloatType> SumProductBP::MarginalOfVar(int WhichVar)
{
    vector< FloatType > RetVal;

    RetVal = EdgeToVarMsg[ Var_Fn_Edge[WhichVar][0] ];
    FloatType Total = 0.0;
    for(int i=0;i<RetVal.size();i++)
    {
        RetVal[i] *= EdgeToFnMsg[ Var_Fn_Edge[WhichVar][0] ][i];
        Total += RetVal[i];
    }
    Total = 1.0 / Total;
    for(int i=0;i<RetVal.size();i++)
    {
        RetVal[i] *= Total;
    }

    return RetVal;
}




void SumProductBP::MarginalOfVars(vector< vector<FloatType> > & Result)
{
    Result.resize(graph.variables.size());
    for(int i=0;i<graph.variables.size();i++)
    {
        Result[i] = EdgeToVarMsg[ Var_Fn_Edge[i][0] ];

        FloatType Total = 0.0;
        for(int j=0;j<Result[i].size();j++)
        {
            Result[i][j] *= EdgeToFnMsg[ Var_Fn_Edge[i][0] ][j];
            Total += Result[i][j];
        }
        Total = 1.0 / Total;
        for(int j=0;j<Result[i].size();j++)
        {
            Result[i][j] *= Total;
        }
    }
}






void SumProductBP::DecisionOfVars(vector< int > & Result)
{
    Result.resize(graph.variables.size());
    for(int i=0;i<graph.variables.size();i++)
    {
        Result[i] = 0 ;

        FloatType MaxProb = -INFINITY;
        for(int j=0;j<graph.variables[i].cardinality;j++)
        {
            FloatType Temp = EdgeToVarMsg[ Var_Fn_Edge[i][0] ][j] * EdgeToFnMsg[ Var_Fn_Edge[i][0] ][j];
            if(Temp > MaxProb)
            {
                MaxProb = Temp;
                Result[i] = j;
            }
        }
    }
}




void SumProductBP::GetDiffDist(vector<FloatType> & theDist)
{
    int MaxD = 0;
    for(auto & vec:EdgeToFnMsg)
    {
        if(vec.size() > MaxD) MaxD = vec.size();
    }
    theDist.assign(MaxD,0.0);
    TempVec.resize(MaxD);
    for(int i=0;i<EdgeToFnMsg.size();i++)
    {
        for(int j=0;j<EdgeToFnMsg[i].size();j++)
        {
            TempVec[j] = abs(EdgeToFnMsg[i][j] - EdgeToFnMsg_old[i][j]);
            // TempVec[j] = TempVec[j]*TempVec[j];
        }
        for(int j=EdgeToFnMsg.size();j<MaxD;j++)
        {
            TempVec[j] = 0.0;
        }
        // sort Descending with reverse iterators
        sort(TempVec.rbegin(),TempVec.rend());
        for(int j=0;j<MaxD;j++)
        {
            theDist[j] += TempVec[j]*TempVec[j];
        }
    }
    // FloatType Total=0.0;
    // for(auto i:theDist)
    //  Total += i;
    // for(auto & i:theDist)
    //  i /= Total;
}





////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////










