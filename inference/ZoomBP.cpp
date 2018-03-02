#include "ZoomBP.h"
#include <cmath>
#include <pthread.h>

#ifndef NULL
#define NULL 0
#endif

int const Levels = 7;
int const DistSize = 10;
FloatType const UpRange = pow( 2.0 * Levels + 1.0 , 20.0 / DistSize );
FloatType const DownRange = pow( 2.0 * Levels + 1.0 , -1.0 / DistSize ) ;
FloatType const RANGE_INIT = 0.0001;
FloatType const RANGE_LB = 1e-19;
FloatType const RANGE_UB = 0.1;

int const ChunkSize = 200;

string ZoomBP::PrintParams()
{
    stringstream ss;
    ss << "Levels=" << 1+2*Levels
        << ",UpRange=" << UpRange
        << ",DownRange=" << DownRange
        << ",RANGE_INIT=" << RANGE_INIT ;
    return ss.str();
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


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// DONE
// Verified by inspection only.
int ZoomBP::UseGraph(FactorGraph newGraph)
{
    GraphSolver::UseGraph(newGraph);

    int EdgeCount = VarOfEdge.size();

    EdgeToFn_NonNormal.resize(EdgeCount);
    EdgeToFn_Dist.resize(EdgeCount);
    EdgeToFn_I.resize(EdgeCount);
    EdgeToFn_D.resize(EdgeCount);
    EdgeToFn_Q.resize(EdgeCount);
    EdgeToFn_R.assign(EdgeCount,RANGE_INIT);
    EdgeToFn_Alpha.resize(EdgeCount);
    EdgeToFn_Gamma.resize(EdgeCount);
    EdgeToFn_DistEst.resize(EdgeCount);
    EdgeToFn_DistEst_old.resize(EdgeCount);

    EdgeToVar_Marginal.resize(EdgeCount);
    EdgeToVar_Dist.resize(EdgeCount);
    EdgeToVar_I.resize(EdgeCount);
    EdgeToVar_D.resize(EdgeCount);
    EdgeToVar_Q.resize(EdgeCount);
    EdgeToVar_R.assign(EdgeCount,RANGE_INIT);
    EdgeToVar_Alpha.resize(EdgeCount);
    EdgeToVar_Gamma.resize(EdgeCount);
    EdgeToVar_DistEst.resize(EdgeCount);
    EdgeToVar_DistEst_old.resize(EdgeCount);

    EdgeToVar_Indices.resize(EdgeCount);

    int WhichEdge = 0;
    // go through each function node
    for(int i=0;i<graph.functions.size();i++)
    {
        // for every variable connected to that function node
        // i.e. every edge from that function node
        for(int j=0;j<graph.functions[i].vars.size();j++)
        {
            int tempSize;
            tempSize = graph.variables[ graph.functions[i].vars[j] ].cardinality;
            // initial assignment of messages is UNIFORM
            EdgeToFn_NonNormal   [WhichEdge].assign( tempSize , 1.0L ) ;
            EdgeToFn_Dist        [WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;
            EdgeToFn_DistEst     [WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;
            EdgeToFn_DistEst_old [WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;
            EdgeToVar_Marginal   [WhichEdge].assign( tempSize , 0.0L ) ;
            EdgeToVar_Dist       [WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;
            EdgeToVar_DistEst    [WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;
            EdgeToVar_DistEst_old[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) ) ;

            EdgeToVar_Indices[WhichEdge].resize(graph.functions[i].vars.size());

            WhichEdge++;
        }
    }


    for(int i=0;i<EdgeToVar_Marginal.size();i++)
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
        for(int j=0;j<EdgeToVar_Marginal[i].size();j++)
        {
            EdgeToVar_Marginal[i][j] = 0.0;
        }

        // zero out the symbols
        for(int j=0;j<EdgeToVar_Indices[i].size();j++)
        {
            EdgeToVar_Indices[i][j] = 0;
        }

        // go through every function value
        for(int j=0;j<graph.functions[theFn].ValuesSize;j++)
        {
            FloatType theProduct = 1.0;
            // product of the incoming messages at the current point (exclude current var)
            for(int k=0;k<EdgeToVar_Indices[i].size();k++)
            {
                // EdgeToFnMsg [desired edge] [value of the variable]
                // desired edge is k'th one attached to the current
                //     function, hence Fn_Var_Edge[theFn][k]
                // value of the variable is a coordinate 'k' from EdgeToVar_Indices[i]
                // note: the vars attached to the fn and the var symbols are same order
                if(k != LeaveOut) theProduct *=
                    EdgeToFn_DistEst[  Fn_Var_Edge[theFn][k]  ][  EdgeToVar_Indices[i][k]  ];
            }
            EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ]
                += graph.values[graph.functions[theFn].values][j] * theProduct;
                // j == CoordsToFullIndex(theFn,EdgeToVar_Indices[i])
            IncSymbols(theFn,EdgeToVar_Indices[i]);
        }

        // normalize the EdgeToVar_Dist
        FloatType MsgTotal = 0.0;
        for(int j=0;j<EdgeToVar_Marginal[i].size();j++)
        {
            MsgTotal += EdgeToVar_Marginal[i][j];
        }
        MsgTotal = 1.0 / MsgTotal;
        EdgeToVar_Dist[i] = EdgeToVar_Marginal[i];
        for(int j=0;j<EdgeToVar_Dist[i].size();j++)
        {
            EdgeToVar_Dist[i][j] *= MsgTotal;
        }
    }

    for(int i=0;i<EdgeToFn_NonNormal.size();i++)
    {
        // "clear" the message to all ones
        for(int k=0;k<EdgeToFn_NonNormal[i].size();k++)
            EdgeToFn_NonNormal[ i ] [ k ] = 1.0;
        // find edge to var messages attached to VarOfEdge[i]
        // go through all the edges/functions connected to the var
        for(int j=0;j<Var_Fn_Edge[ VarOfEdge[i] ].size();j++)
        {
            // one of them should be the current edge, don't use in update
            if(Var_Fn_Edge[ VarOfEdge[i] ][j] != i)
                for(int k=0;k<EdgeToFn_NonNormal[i].size();k++)
                    EdgeToFn_NonNormal[ i ] [k] *=
                        EdgeToVar_DistEst[ Var_Fn_Edge[ VarOfEdge[i] ][j] ] [k];
        }

        // normalize the EdgeToFn_Dist
        FloatType MsgTotal = 0.0;
        for(int j=0;j<EdgeToFn_NonNormal[i].size();j++)
        {
            MsgTotal += EdgeToFn_NonNormal[i][j];
        }
        MsgTotal = 1.0 / MsgTotal;
        EdgeToFn_Dist[i] = EdgeToFn_NonNormal[i];
        for(int j=0;j<EdgeToFn_Dist[i].size();j++)
        {
            EdgeToFn_Dist[i][j] *= MsgTotal;
        }
    }

    return 0;
}




















void ZoomBP::ResetState()
{
    // EdgeToFn_DistEst - uniform
    // EdgeToFn_Dist - uniform
    // EdgeToVar_Marginal - must be computed
    // EdgeToVar_DistEst - uniform
    // EdgeToFn_R
    // EdgeToVar_R
    for(int i=0;i<EdgeToFn_DistEst.size();i++)
    {
        int TempSize = EdgeToFn_DistEst[i].size();
        for(int j=0;j<TempSize;j++)
        {
            EdgeToFn_DistEst[i][j] =
            EdgeToFn_Dist[i][j] =
            EdgeToVar_DistEst[i][j] =
                1.0 / (FloatType)TempSize;
        }
        EdgeToFn_R[i] = EdgeToVar_R[i] = RANGE_INIT;
    }

    for(int i=0;i<EdgeToVar_Marginal.size();i++)
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
        for(int j=0;j<EdgeToVar_Marginal[i].size();j++)
        {
            EdgeToVar_Marginal[i][j] = 0.0;
        }

        // zero out the symbols
        for(int j=0;j<EdgeToVar_Indices[i].size();j++)
        {
            EdgeToVar_Indices[i][j] = 0;
        }

        // go through every function value
        for(int j=0;j<graph.functions[theFn].ValuesSize;j++)
        {
            FloatType theProduct = 1.0;
            // product of the incoming messages at the current point (exclude current var)
            for(int k=0;k<EdgeToVar_Indices[i].size();k++)
            {
                // EdgeToFnMsg [desired edge] [value of the variable]
                // desired edge is k'th one attached to the current
                //     function, hence Fn_Var_Edge[theFn][k]
                // value of the variable is a coordinate 'k' from EdgeToVar_Indices[i]
                // note: the vars attached to the fn and the var symbols are same order
                if(k != LeaveOut) theProduct *=
                    EdgeToFn_DistEst[  Fn_Var_Edge[theFn][k]  ][  EdgeToVar_Indices[i][k]  ];
            }
            EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ]
                += graph.values[graph.functions[theFn].values][j] * theProduct;
                // j == CoordsToFullIndex(theFn,EdgeToVar_Indices[i])
            IncSymbols(theFn,EdgeToVar_Indices[i]);
        }
    }
}


void ZoomBP::GetState(vector< vector< FloatType > > & State)
{
    // State consists of:
    //  - EdgeToFn_Dist
    //  - EdgeToFn_DistEst
    //  - EdgeToFn_R
    //  - EdgeToVar_DistEst
    //  - EdgeToVar_R
    State.resize(EdgeToFn_Dist.size()+EdgeToFn_DistEst.size()+1+EdgeToVar_DistEst.size()+1);
    int Index=0;
    for(auto & v:EdgeToFn_Dist)
    {
        State[Index++] = v;
    }
    for(auto & v:EdgeToFn_DistEst)
    {
        State[Index++] = v;
    }
    State[Index++] = EdgeToFn_R;
    for(auto & v:EdgeToVar_DistEst)
    {
        State[Index++] = v;
    }
    State[Index++] = EdgeToVar_R;
}


void ZoomBP::RestoreState(vector< vector< FloatType > > const & State)
{
    int Index=0;
    for(auto & v:EdgeToFn_Dist)
    {
        v = State[Index++];
    }
    for(auto & v:EdgeToFn_DistEst)
    {
        v = State[Index++];
    }
    EdgeToFn_R = State[Index++];
    for(auto & v:EdgeToVar_DistEst)
    {
        v = State[Index++];
    }
    EdgeToVar_R = State[Index++];
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






/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void ZoomBP::Iter()
{
    // update every EdgeToFnMsg
    for(int i=0;i<EdgeToFn_Dist.size();i++)
    {
        UpdateMsgFromVar_Edge(i);
    }

    // update every EdgeToVarMsg
    for(int i=0;i<EdgeToVar_Dist.size();i++)
    {
        UpdateMsgFromFn_Edge(i);
    }

    // update every EdgeToVarMsg
    for(int i=0;i<EdgeToFn_Dist.size();i++)
    {
        UpdateBeliefFromVar_Edge(i);
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


///////////////////////
struct ZoomBP::ThreadData
{
    ZoomBP * Graph;
    int Progress;
    pthread_mutex_t Mutex;
};
///////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void * ZoomBP::EdgeToFnThread(void * d)
{
    // cast the input pointer so we can use it
    ThreadData * theData = (ThreadData *)d;

    int QuitThread;
    do
    {
        int StartCompute;
        QuitThread = 0;
        pthread_mutex_lock(&(theData->Mutex));
        if(theData->Progress < theData->Graph->EdgeToFn_Dist.size())
        {
            StartCompute = theData->Progress;
            theData->Progress += ChunkSize;
        }
        else QuitThread = 1;
        pthread_mutex_unlock(&(theData->Mutex));
        if(!QuitThread)
        {
            for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToFn_Dist.size();i++)
                theData->Graph->UpdateMsgFromVar_Edge(i);
        }
    }while(!QuitThread);
    return NULL;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void * ZoomBP::EdgeToVarThread(void * d)
{
    // cast the input pointer so we can use it
    ThreadData * theData = (ThreadData *)d;

    int QuitThread;
    do
    {
        int StartCompute;
        QuitThread = 0;
        pthread_mutex_lock(&(theData->Mutex));
        if(theData->Progress < theData->Graph->EdgeToVar_Dist.size())
        {
            StartCompute = theData->Progress;
            theData->Progress += ChunkSize;
        }
        else QuitThread = 1;
        pthread_mutex_unlock(&(theData->Mutex));
        if(!QuitThread)
        {
            for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToVar_Dist.size();i++)
                theData->Graph->UpdateMsgFromFn_Edge(i);
        }
    }while(!QuitThread);
    return NULL;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void * ZoomBP::BeliefFromVarThread(void * d)
{
    // cast the input pointer so we can use it
    ThreadData * theData = (ThreadData *)d;

    int QuitThread;
    do
    {
        int StartCompute;
        QuitThread = 0;
        pthread_mutex_lock(&(theData->Mutex));
        if(theData->Progress < theData->Graph->EdgeToFn_Dist.size())
        {
            StartCompute = theData->Progress;
            theData->Progress += ChunkSize;
        }
        else QuitThread = 1;
        pthread_mutex_unlock(&(theData->Mutex));
        if(!QuitThread)
        {
            for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToFn_Dist.size();i++)
                theData->Graph->UpdateBeliefFromVar_Edge(i);
        }
    }while(!QuitThread);
    return NULL;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void ZoomBP::IterThreads(int NumThreads)
{
    if(NumThreads==1)
    {
        Iter();
        return;
    }
    //////////////////////////////////////////////
    // update every EdgeToFn Dist
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
    // update every EdgeToVar Dist
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
    // update every EdgeToFn Belief
    {
        vector<pthread_t> myThreads;
        myThreads.resize(NumThreads);

        ThreadData theData;

        theData.Graph = this;
        theData.Progress = 0;
        pthread_mutex_init(&theData.Mutex,NULL);

        for(int i=0;i<NumThreads;i++)
        {
            pthread_create( &myThreads[i], NULL, BeliefFromVarThread, (void*)(&theData) );
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






/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/*
 * update the (optionally quantized) var to fn message
 *   - update index
 *   - quantized step value (corresponding with R), or step size (direction is implicit)
 *
 * update alpha, gamma, and the estimate of the var to fn distribution
 *
 */
void ZoomBP::UpdateMsgFromVar_Edge(int i)
{
    EdgeToFn_DistEst_old[i] = EdgeToFn_DistEst[i];

    // the "message"
    EdgeToFn_I[i] = 0;

    ///////////////////////////////////////////////////////
    // determine 'EdgeToFn_I[i]' - which entry to update
    ///////////////////////////////////////////////////////
    FloatType Prod1 = 0.0;
    FloatType MagSq = 1.0;
    for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
    {
        Prod1 += EdgeToFn_DistEst[i][j] * (EdgeToFn_DistEst[i][j] - EdgeToFn_Dist[i][j]);
        MagSq += EdgeToFn_DistEst[i][j] * EdgeToFn_DistEst[i][j];
    }
    FloatType MaxMetric = 0.0;
    for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
    {
        FloatType tempMetric = Prod1 + EdgeToFn_Dist[i][j] - EdgeToFn_DistEst[i][j];
        tempMetric = abs(tempMetric) / (sqrt(MagSq - 2.0 * EdgeToFn_DistEst[i][j]));
        if( MaxMetric < tempMetric )
        {
            MaxMetric = tempMetric;
            EdgeToFn_I[i] = j;
        }
    }

    //////////////////////////////////////////////////////////
    // determine 'EdgeToFn_Q[i]' - which quantization value
    //////////////////////////////////////////////////////////
    FloatType OldEst = EdgeToFn_DistEst[i][EdgeToFn_I[i]];
    FloatType Projection = Prod1 + EdgeToFn_Dist[i][EdgeToFn_I[i]] - OldEst;
    Projection = Projection / (sqrt(MagSq - 2.0 * OldEst));
    if(QuantVarToFn)
    {
        Projection *= (2.0*Levels+1.0) / EdgeToFn_R[i];
        EdgeToFn_Q[i] = (int)(round(Projection));
        if(EdgeToFn_Q[i] < -Levels) EdgeToFn_Q[i] = -Levels;
        if(EdgeToFn_Q[i] >  Levels) EdgeToFn_Q[i] =  Levels;
    }

    /////////////////////////////////////////////////////////////////////
    // send EdgeToFn_I[i] and EdgeToFn_Q[i] to receiving node......
    /////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////
    // compute estimate update
    //////////////////////////////////////////////////////
    if(QuantVarToFn)
    {
        EdgeToFn_DistEst[i][EdgeToFn_I[i]] +=
            (FloatType)(EdgeToFn_Q[i]) * EdgeToFn_R[i] / (2.0*Levels+1.0)
          / (sqrt(MagSq - 2.0 * OldEst))
          * ( 1.0 - OldEst );
    }
    else
    {
        EdgeToFn_DistEst[i][EdgeToFn_I[i]] +=
            Projection
          / (sqrt(MagSq - 2.0 * OldEst))
          * ( 1.0 - OldEst );
    }
    if(EdgeToFn_DistEst[i][EdgeToFn_I[i]] <= 0.0) EdgeToFn_DistEst[i][EdgeToFn_I[i]] = OldEst/(2.0*Levels+1.0);
    if(EdgeToFn_DistEst[i][EdgeToFn_I[i]] >= 1.0) EdgeToFn_DistEst[i][EdgeToFn_I[i]] = 1.0 - (1.0-OldEst)/(2.0*Levels+1.0);

    FloatType Delta = EdgeToFn_DistEst[i][EdgeToFn_I[i]] - OldEst;
    EdgeToFn_D[i] = Delta;

    FloatType Alpha = 1.0 - Delta / (1.0 - OldEst ) ;
    EdgeToFn_Alpha[i] = Alpha;

    FloatType Gamma = ( 1.0 - Alpha ) * OldEst + Delta ;
    EdgeToFn_Gamma[i] = Gamma;

    FloatType Total = 0.0;
    for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
    {
        if(j != EdgeToFn_I[i]) EdgeToFn_DistEst[i][j] *= Alpha;
        Total += EdgeToFn_DistEst[i][j];
    }
    for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
    {
        EdgeToFn_DistEst[i][j] *= (1.0/Total);
    }

    // compute state update
    if(QuantVarToFn)
    {
        if( EdgeToFn_Q[i] == -Levels || EdgeToFn_Q[i] == Levels ) EdgeToFn_R[i] *= UpRange;
        else EdgeToFn_R[i] *= DownRange;
        if(EdgeToFn_R[i] < RANGE_LB) EdgeToFn_R[i] = RANGE_LB;
        if(EdgeToFn_R[i] > RANGE_UB) EdgeToFn_R[i] = RANGE_UB;
    }
}




/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


/*
 * update the fn to var output marginal distribution
 *
 * update the (optionally quantized) fn to var message
 *   - update index
 *   - quantized step value (corresponding with R), or step size (direction is implicit)
 *
 * update alpha, gamma, and the estimate of the fn to var distribution
 *
 */
void ZoomBP::UpdateMsgFromFn_Edge(int i)
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


    // update EdgeToVar_Marginal[i]
    // this is done with a series of updates of the form:
    //  dist(l) = alpha * dist(l) + gamma * sum_(not leaveout or modindex)[ A(...)*prod[dist_est] ]
    // one update per incoming other edge (degree-1) i.e. exclude LeaveOut
    for(int j=0;j<graph.functions[theFn].vars.size();j++)
    {
        if(j!=LeaveOut)
        {
            // update for message from variable 'j'
            // this uses edge number Fn_Var_Edge[theFn][j]

            // zero out the symbols
            for(int k=0;k<EdgeToVar_Indices[i].size();k++)
            {
                EdgeToVar_Indices[i][k] = 0;
            }
            // summation index for the updating distribution should always be the same
            //  * we are now incorporating the update for distribution 'j'
            EdgeToVar_Indices[i][j] = EdgeToFn_I[Fn_Var_Edge[theFn][j]];
            // Scale down the marginal by Alpha
            for(int k=0;k<EdgeToVar_Marginal[i].size();k++)
            {
                EdgeToVar_Marginal[i][k] *= EdgeToFn_Alpha[Fn_Var_Edge[theFn][j]] ;
            }

            for(int k=0;k<graph.functions[theFn].ValuesSize / graph.variables[graph.functions[theFn].vars[j]].cardinality ; k++)
            {
                // this product should leave out 'LeaveOut' as well as 'j'
                // 'LeaveOut' is certainly excluded
                // 'j' is excluded because it is accounted for later in Gamma.
                // vars below 'j' should take new dist estimate
                // vars after 'j' should take old dist estimate
                FloatType ProductOfInputs = 1.0 ;
                for(int l=0;l<graph.functions[theFn].vars.size();l++)
                {
                    if(l != LeaveOut)
                    {
                        // use new dist estimate
                        if(l < j) ProductOfInputs *= EdgeToFn_DistEst[ i ][ EdgeToVar_Indices[i][l] ];
                        // use old dist estimate
                        if(l > j) ProductOfInputs *= EdgeToFn_DistEst_old[ i ][ EdgeToVar_Indices[i][l] ];
                    }
                }
                EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ] +=
                    graph.values[graph.functions[theFn].values][ CoordsToFullIndex(theFn,EdgeToVar_Indices[i]) ]
                    * ProductOfInputs
                    * EdgeToFn_Gamma[Fn_Var_Edge[theFn][j]] ;

                IncSymbolsExcept(theFn,EdgeToVar_Indices[i],j);
            }
        }
    }

    // update EdgeToVar_Dist[i]
    // just a normalized EdgeToVar_Marginal[i]
    {
        FloatType Total = 0.0;
        for(int j=0;j<EdgeToVar_Marginal[i].size();j++)
        {
            Total += EdgeToVar_Marginal[i][j];
        }
        Total = 1.0 / Total;
        for(int j=0;j<EdgeToVar_Marginal[i].size();j++)
        {
            EdgeToVar_Dist[i][j] = EdgeToVar_Marginal[i][j] * Total;
        }
    }

    {
        EdgeToVar_DistEst_old[i] = EdgeToVar_DistEst[i];

        // the "message"
        EdgeToVar_I[i] = 0;

        ////////////////////////////////////////////////////////
        // determine 'EdgeToVar_I[i]' - which entry to update
        ////////////////////////////////////////////////////////
        FloatType Prod1 = 0.0;
        FloatType MagSq = 1.0;
        for( int j=0 ; j<EdgeToVar_DistEst[i].size() ; j++ )
        {
            Prod1 += EdgeToVar_DistEst[i][j] * (EdgeToVar_DistEst[i][j] - EdgeToVar_Dist[i][j]);
            MagSq += EdgeToVar_DistEst[i][j] * EdgeToVar_DistEst[i][j];
        }
        FloatType MaxMetric = 0.0;
        for( int j=0 ; j<EdgeToVar_DistEst[i].size() ; j++ )
        {
            FloatType tempMetric = Prod1 + EdgeToVar_Dist[i][j] - EdgeToVar_DistEst[i][j];
            tempMetric = abs(tempMetric) / (sqrt(MagSq - 2.0 * EdgeToVar_DistEst[i][j]));
            if( MaxMetric < tempMetric )
            {
                MaxMetric = tempMetric;
                EdgeToVar_I[i] = j;
            }
        }

        ///////////////////////////////////////////////////////////
        // determine 'EdgeToVar_Q[i]' - which quantization value
        ///////////////////////////////////////////////////////////
        FloatType OldEst = EdgeToVar_DistEst[i][EdgeToVar_I[i]];
        FloatType Projection = Prod1 + EdgeToVar_Dist[i][EdgeToVar_I[i]] - OldEst;
        Projection /= (sqrt(MagSq - 2.0 * OldEst));
        if(QuantFnToVar)
        {
            Projection *= (2.0*Levels+1.0) / EdgeToVar_R[i];
            EdgeToVar_Q[i] = (int)(round(Projection));
            if(EdgeToVar_Q[i] < -Levels) EdgeToVar_Q[i] = -Levels;
            if(EdgeToVar_Q[i] >  Levels) EdgeToVar_Q[i] =  Levels;
        }

        ///////////////////////////////////////////////////////////////////////
        // send EdgeToVar_I[i] and EdgeToVar_Q[i] to receiving node......
        ///////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////
        // compute estimate update
        //////////////////////////////////////////////////////
        if(QuantFnToVar)
        {
            EdgeToVar_DistEst[i][EdgeToVar_I[i]] +=
                (FloatType)(EdgeToVar_Q[i]) * EdgeToVar_R[i] / (2.0*Levels+1.0)
              / (sqrt(MagSq - 2.0 * OldEst))
              * ( 1.0 - OldEst );
        }
        else
        {
            EdgeToVar_DistEst[i][EdgeToVar_I[i]] +=
                Projection
              / (sqrt(MagSq - 2.0 * OldEst))
              * ( 1.0 - OldEst );
        }
        if(EdgeToVar_DistEst[i][EdgeToVar_I[i]] <= 0.0) EdgeToVar_DistEst[i][EdgeToVar_I[i]] = OldEst/(2.0*Levels+1.0);
        if(EdgeToVar_DistEst[i][EdgeToVar_I[i]] >= 1.0) EdgeToVar_DistEst[i][EdgeToVar_I[i]] = 1.0 - (1.0-OldEst)/(2.0*Levels+1.0);

        FloatType Delta = EdgeToVar_DistEst[i][EdgeToVar_I[i]] - OldEst;
        EdgeToVar_D[i] = Delta;

        FloatType Alpha = 1.0 - Delta / (1.0 - OldEst ) ;
        EdgeToVar_Alpha[i] = Alpha;

        FloatType Gamma = ( 1.0 - Alpha ) * OldEst + Delta ;
        EdgeToVar_Gamma[i] = Gamma;

        FloatType Total = 0.0;
        for(int j=0;j<EdgeToVar_DistEst[i].size();j++)
        {
            if(j != EdgeToVar_I[i]) EdgeToVar_DistEst[i][j] *= Alpha;
            Total += EdgeToVar_DistEst[i][j];
        }
        for(int j=0;j<EdgeToVar_DistEst[i].size();j++)
        {
            EdgeToVar_DistEst[i][j] *= (1.0/Total);
        }

        // compute state update
        if(QuantFnToVar)
        {
            if( EdgeToVar_Q[i] == -Levels || EdgeToVar_Q[i] == Levels ) EdgeToVar_R[i] *= UpRange;
            else EdgeToVar_R[i] *= DownRange;
            if(EdgeToVar_R[i] < RANGE_LB) EdgeToVar_R[i] = RANGE_LB;
            if(EdgeToVar_R[i] > RANGE_UB) EdgeToVar_R[i] = RANGE_UB;
        }
    }
}





/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/*
 * update the var to fn output belief distribution
 *
 */
void ZoomBP::UpdateBeliefFromVar_Edge(int i)
{
    // ALTERNATIVELY, THIS COULD BE A LOWER COMPLEXITY UPDATE

    // "clear" the message to all ones
    for(int k=0;k<EdgeToFn_NonNormal[ i ].size();k++)
        EdgeToFn_NonNormal[ i ] [ k ] = 1.0;
    // find edge to var messages attached to VarOfEdge[i]
    // go through all the edges/functions connected to the var
    for(int j=0;j<Var_Fn_Edge[ VarOfEdge[i] ].size();j++)
    {
        // one of them should be the current edge, don't use in update
        if(Var_Fn_Edge[ VarOfEdge[i] ][j] != i)
            for(int k=0;k<EdgeToFn_NonNormal[ i ].size();k++)
                EdgeToFn_NonNormal[ i ] [ k ] *=
                    EdgeToVar_DistEst[ Var_Fn_Edge[ VarOfEdge[i] ][j] ] [k];
    }

    // normalize the EdgeToFn_Dist
    FloatType MsgTotal = 0.0;
    for(int j=0;j<EdgeToFn_NonNormal[i].size();j++)
    {
        MsgTotal += EdgeToFn_NonNormal[i][j];
    }
    MsgTotal = 1.0 / MsgTotal;
    EdgeToFn_Dist[i] = EdgeToFn_NonNormal[i];
    for(int j=0;j<EdgeToFn_Dist[i].size();j++)
    {
        EdgeToFn_Dist[i][j] *= MsgTotal;
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





/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


// DONE
vector<FloatType> ZoomBP::MarginalOfVar(int WhichVar)
{
    vector< FloatType > RetVal;

    RetVal = EdgeToVar_Dist[ Var_Fn_Edge[WhichVar][0] ];
    FloatType Total = 0.0;
    for(int i=0;i<RetVal.size();i++)
    {
        RetVal[i] *= EdgeToFn_Dist[ Var_Fn_Edge[WhichVar][0] ][i];
        Total += RetVal[i];
    }
    Total = 1.0 / Total;
    for(int i=0;i<RetVal.size();i++)
    {
        RetVal[i] *= Total;
    }

    return RetVal;
}




void ZoomBP::MarginalOfVars(vector< vector<FloatType> > & Result)
{
    Result.resize(graph.variables.size());
    for(int i=0;i<graph.variables.size();i++)
    {
        Result[i] = EdgeToVar_Dist[ Var_Fn_Edge[i][0] ];

        FloatType Total = 0.0;
        for(int j=0;j<Result[i].size();j++)
        {
            Result[i][j] *= EdgeToFn_Dist[ Var_Fn_Edge[i][0] ][j];
            Total += Result[i][j];
        }
        Total = 1.0 / Total;
        for(int j=0;j<Result[i].size();j++)
        {
            Result[i][j] *= Total;
        }
    }
}





void ZoomBP::DecisionOfVars(vector< int > & Result)
{
    Result.resize(graph.variables.size());
    for(int i=0;i<graph.variables.size();i++)
    {
        Result[i] = 0 ;

        FloatType MaxProb = -INFINITY;
        for(int j=0;j<graph.variables[i].cardinality;j++)
        {
            FloatType Temp = EdgeToVar_Dist[ Var_Fn_Edge[i][0] ][j] * EdgeToFn_Dist[ Var_Fn_Edge[i][0] ][j];
            if(Temp > MaxProb)
            {
                MaxProb = Temp;
                Result[i] = j;
            }
        }
    }
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










