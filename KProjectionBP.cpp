#include "KProjectionBP.h"
#include <cmath>
#include <pthread.h>
#include <cstdlib>
#include <cassert>

#ifndef NULL
#define NULL 0
#endif

#define NOCHECKS

int const ChunkSize = 200;

string KProjectionBP::PrintParams()
{
	stringstream ss;
	ss << "KUpdate=" << KUpdate ;
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
int KProjectionBP::UseGraph(FactorGraph newGraph)
{
	GraphSolver::UseGraph(newGraph);

	int EdgeCount = VarOfEdge.size();

	EdgeToFn_NonNormal.resize(EdgeCount);
	EdgeToFn_Dist.resize(EdgeCount);
	EdgeToFn_Alpha.resize(EdgeCount);
	EdgeToFn_Gamma.resize(EdgeCount);
	EdgeToFn_DistEst.resize(EdgeCount);
	EdgeToFn_DistEst_old.resize(EdgeCount);
	EdgeToFn_FindK.resize(EdgeCount);

	EdgeToVar_Marginal.resize(EdgeCount);
	EdgeToVar_Dist.resize(EdgeCount);
	EdgeToVar_Alpha.resize(EdgeCount);
	EdgeToVar_Gamma.resize(EdgeCount);
	EdgeToVar_DistEst.resize(EdgeCount);
	EdgeToVar_DistEst_old.resize(EdgeCount);
	EdgeToVar_FindK.resize(EdgeCount);

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

			EdgeToFn_FindK[WhichEdge].resize(tempSize);
			EdgeToVar_FindK[WhichEdge].resize(tempSize);
			for(int k=0;k<tempSize;k++)
			{
				EdgeToFn_FindK[WhichEdge][k] = EdgeToVar_FindK[WhichEdge][k] = k;
			}

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












void KProjectionBP::ResetState()
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



void KProjectionBP::GetState(vector< vector< FloatType > > & State)
{
	// State consists of:
	//  - EdgeToFn_Dist[][] 
	//  - EdgeToFn_DistEst[][]
	//  - EdgeToVar_Marginal[][] (it could be recomputed from EdgeToFn_DistEst...)
	//  - EdgeToVar_DistEst[][]
	//  - EdgeToVar_DistEst_old[][]
	State.resize(EdgeToFn_Dist.size()+EdgeToFn_DistEst.size()+EdgeToVar_Marginal.size()+EdgeToVar_DistEst.size()+EdgeToVar_DistEst_old.size());
	int Index=0;
	for(auto & v:EdgeToFn_Dist)
	{
		State[Index++] = v;
	}
	for(auto & v:EdgeToFn_DistEst)
	{
		State[Index++] = v;
	}
	for(auto & v:EdgeToVar_Marginal)
	{
		State[Index++] = v;
	}
	for(auto & v:EdgeToVar_DistEst)
	{
		State[Index++] = v;
	}
	for(auto & v:EdgeToVar_DistEst_old)
	{
		State[Index++] = v;
	}
}


void KProjectionBP::RestoreState(vector< vector< FloatType > > const & State)
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
	for(auto & v:EdgeToVar_Marginal)
	{
		v = State[Index++];
	}
	for(auto & v:EdgeToVar_DistEst)
	{
		v = State[Index++];
	}
	for(auto & v:EdgeToVar_DistEst_old)
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

void KProjectionBP::Iter()
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
struct KProjectionBP::ThreadData
{
	KProjectionBP * Graph;
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

void * KProjectionBP::EdgeToFnThread(void * d)
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


void * KProjectionBP::EdgeToVarThread(void * d)
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


void * KProjectionBP::BeliefFromVarThread(void * d)
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


void KProjectionBP::IterThreads(int NumThreads)
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

FloatType const MinThreshold = 1e-15;
void KProjectionBP::UpdateMsgFromVar_Edge(int i)
{
	///////////////////////////////////////////////////////////
	// determine 'EdgeToFn_FindK[i]' - which entries to update
	///////////////////////////////////////////////////////////
	for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
	{
		// temporary use of EdgeToFn_DistEst_old[i], rather than allocate another
		EdgeToFn_DistEst_old[i][j] = abs( EdgeToFn_Dist[i][j] - EdgeToFn_DistEst[i][j] );
	}
	// the update indices will be in positions 0 thru KUpdate-1 in EdgeToFn_FindK[i]
	quick_select_largest(EdgeToFn_DistEst_old[i],EdgeToFn_FindK[i],0,EdgeToFn_DistEst[i].size()-1,KUpdate);

	//////////////////////////////////////////////////////////////
	// now use EdgeToFn_DistEst_old[i] for its intended purpose
	//////////////////////////////////////////////////////////////
	EdgeToFn_DistEst_old[i] = EdgeToFn_DistEst[i];

	//////////////////////////////////////////////////////
	// compute estimate update
	//////////////////////////////////////////////////////
	FloatType DeltaSum = 0.0;
	FloatType OldEstSum = 0.0;
	for(int j=0;j<KUpdate;j++)
	{
		// estimate equals actual value in select positions
		EdgeToFn_DistEst[ i ][ EdgeToFn_FindK[i][j] ] = EdgeToFn_Dist[ i ][ EdgeToFn_FindK[i][j] ];
		DeltaSum += EdgeToFn_DistEst[ i ][ EdgeToFn_FindK[i][j] ] - EdgeToFn_DistEst_old[ i ][ EdgeToFn_FindK[i][j] ];
		OldEstSum += EdgeToFn_DistEst_old[ i ][ EdgeToFn_FindK[i][j] ] ;
	}

	FloatType MinEst = 1.0;
	for(int j=KUpdate;j<EdgeToFn_DistEst[i].size();j++)
	{
		if(MinEst > EdgeToFn_DistEst[i][ EdgeToFn_FindK[i][j] ])
			MinEst = EdgeToFn_DistEst[i][ EdgeToFn_FindK[i][j] ];
	}

	FloatType Alpha = 1.0 - DeltaSum / (1.0 - OldEstSum ) ;
	EdgeToFn_Alpha[i] = Alpha;
	if(Alpha*MinEst > MinThreshold)
	{
		// scale by alpha, also keep track of total
		FloatType Total = 0.0;
		for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
		{
			if(j>=KUpdate)
				EdgeToFn_DistEst[i][ EdgeToFn_FindK[i][j] ] *= Alpha;
			Total += EdgeToFn_DistEst[i][ EdgeToFn_FindK[i][j] ];
		}

		// renormalize for good measure (shouldn't be necessary, could try to take it out)
		for(int j=0;j<EdgeToFn_DistEst[i].size();j++)
		{
			EdgeToFn_DistEst[i][j] *= (1.0/Total);
		}
	}
	else
	{
		EdgeToFn_DistEst[i] = EdgeToFn_DistEst_old[i];
		Alpha = 1.0;
		EdgeToFn_Alpha[i] = Alpha;
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
void KProjectionBP::UpdateMsgFromFn_Edge(int i)
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

#ifndef NOCHECKS
		{
			vector< FloatType > const & Numbers = EdgeToVar_Marginal[i];
			for(int i=0;i<Numbers.size();i++)
				assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
		}
#endif

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
			int WhichIncomingEdge = Fn_Var_Edge[theFn][j];

			if(EdgeToFn_Alpha[WhichIncomingEdge] == 1.0) continue;

#ifndef NOCHECKS
			{
				vector< FloatType > const & Numbers = EdgeToVar_Marginal[i];
				for(int i=0;i<Numbers.size();i++)
				{
					if(isnan(Numbers[i]) || isinf(Numbers[i])) cerr << Numbers[i] << endl;
					assert(!isnan(Numbers[i]));
					assert(!isinf(Numbers[i]));
				}
			}
#endif

			// Scale down the marginal by Alpha
			for(int k=0;k<EdgeToVar_Marginal[i].size();k++)
			{
				EdgeToVar_Marginal[i][k] *= EdgeToFn_Alpha[WhichIncomingEdge] ;
			}

#ifndef NOCHECKS
			{
				vector< FloatType > const & Numbers = EdgeToVar_Marginal[i];
				for(int i=0;i<Numbers.size();i++)
				{
					if(isnan(Numbers[i]) || isinf(Numbers[i])) cerr << Numbers[i] << endl;
					assert(!isnan(Numbers[i]));
					assert(!isinf(Numbers[i]));
				}
			}
#endif

			// zero out the symbols
			for(int k=0;k<EdgeToVar_Indices[i].size();k++)
			{
				EdgeToVar_Indices[i][k] = 0;
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
					if(l == LeaveOut) continue;

					// use new dist estimate
					if(l < j) ProductOfInputs *= EdgeToFn_DistEst[ i ][ EdgeToVar_Indices[i][l] ];
					// use old dist estimate
					if(l > j) ProductOfInputs *= EdgeToFn_DistEst_old[ i ][ EdgeToVar_Indices[i][l] ];
				}
				for( int WhichUpdate = 0 ; WhichUpdate<KUpdate ; WhichUpdate++ )
				{
					// summation index for the updating distribution should always be the same
					//  * we are now incorporating the update for distribution 'j'
					EdgeToVar_Indices[i][j] = EdgeToFn_FindK[WhichIncomingEdge][WhichUpdate];
					// new - alpha old
					FloatType Gamma =
						EdgeToFn_DistEst[WhichIncomingEdge][EdgeToFn_FindK[WhichIncomingEdge][WhichUpdate]]
						-	EdgeToFn_Alpha[WhichIncomingEdge]
							* EdgeToFn_DistEst_old[WhichIncomingEdge][EdgeToFn_FindK[WhichIncomingEdge][WhichUpdate]];

					EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ] +=
						graph.values[graph.functions[theFn].values][ CoordsToFullIndex(theFn,EdgeToVar_Indices[i]) ]
						* ProductOfInputs
						* Gamma ;

#ifndef NOCHECKS
					assert(!isnan(graph.values[graph.functions[theFn].values][ CoordsToFullIndex(theFn,EdgeToVar_Indices[i]) ]));
					assert(!isnan(ProductOfInputs));
					assert(!isnan(Gamma));
					if(isnan(EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ]))
					{
						cerr << graph.values[graph.functions[theFn].values][ CoordsToFullIndex(theFn,EdgeToVar_Indices[i]) ] << endl;
						cerr << ProductOfInputs << endl;
						cerr << Gamma << endl;
						cerr << EdgeToFn_DistEst[WhichIncomingEdge][EdgeToFn_FindK[WhichIncomingEdge][WhichUpdate]] << endl;
						cerr << EdgeToFn_Alpha[WhichIncomingEdge] << endl;
						cerr << EdgeToFn_DistEst_old[WhichIncomingEdge][EdgeToFn_FindK[WhichIncomingEdge][WhichUpdate]] << endl;
					}
					assert(!isnan(EdgeToVar_Marginal[i][ EdgeToVar_Indices[i][LeaveOut] ]));
#endif
				}

				IncSymbolsExcept(theFn,EdgeToVar_Indices[i],j);
			}
		}
	}

	// update EdgeToVar_Dist[i]
	// just a normalized EdgeToVar_Marginal[i]
	{
#ifndef NOCHECKS
		{
			vector< FloatType > const & Numbers = EdgeToVar_Marginal[i];
			for(int i=0;i<Numbers.size();i++)
				assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
		}
#endif

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
		///////////////////////////////////////////////////////
		// determine 'EdgeToFn_I[i]' - which entry to update
		///////////////////////////////////////////////////////
		for(int j=0;j<EdgeToVar_DistEst[i].size();j++)
		{
			// temporary use of EdgeToFn_DistEst_old[i], rather than allocate another
			EdgeToVar_DistEst_old[i][j] = abs( EdgeToVar_Dist[i][j] - EdgeToVar_DistEst[i][j] );
		}
		// the update indices will be in positions 0 thru KUpdate-1 in EdgeToFn_FindK[i]
#ifndef NOCHECKS
		{
			vector< FloatType > const & Numbers = EdgeToVar_DistEst_old[i];
			for(int i=0;i<Numbers.size();i++)
				assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
		}
#endif
		quick_select_largest(EdgeToVar_DistEst_old[i],EdgeToVar_FindK[i],0,EdgeToVar_DistEst[i].size()-1,KUpdate);

		//////////////////////////////////////////////////////////////
		// now use EdgeToFn_DistEst_old[i] for its intended purpose
		//////////////////////////////////////////////////////////////
		EdgeToVar_DistEst_old[i] = EdgeToVar_DistEst[i];

		//////////////////////////////////////////////////////
		// compute estimate update
		//////////////////////////////////////////////////////
		FloatType DeltaSum = 0.0;
		FloatType OldEstSum = 0.0;
#ifndef NOCHECKS
		{
			vector< FloatType > const & Numbers = EdgeToVar_DistEst[i];
			for(int i=0;i<Numbers.size();i++)
				assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
		}
#endif

		for(int j=0;j<KUpdate;j++)
		{
			// estimate equals actual value in select positions
			EdgeToVar_DistEst[ i ][ EdgeToVar_FindK[i][j] ] = EdgeToVar_Dist[ i ][ EdgeToVar_FindK[i][j] ];
			DeltaSum += EdgeToVar_DistEst[ i ][ EdgeToVar_FindK[i][j] ] - EdgeToVar_DistEst_old[ i ][ EdgeToVar_FindK[i][j] ];
			OldEstSum += EdgeToVar_DistEst_old[ i ][ EdgeToVar_FindK[i][j] ] ;
		}

		FloatType MinEst = 1.0;
		for(int j=KUpdate;j<EdgeToVar_DistEst[i].size();j++)
		{
			if(MinEst > EdgeToVar_DistEst[i][ EdgeToVar_FindK[i][j] ])
				MinEst = EdgeToVar_DistEst[i][ EdgeToVar_FindK[i][j] ];
		}

		FloatType Alpha = 1.0 - DeltaSum / (1.0 - OldEstSum ) ;
		EdgeToVar_Alpha[i] = Alpha;
		if(Alpha*MinEst > MinThreshold)
		{
			// scale by alpha, also keep track of total
			FloatType Total = 0.0;
			for(int j=0;j<EdgeToVar_DistEst[i].size();j++)
			{
				if(j>=KUpdate)
					EdgeToVar_DistEst[i][ EdgeToVar_FindK[i][j] ] *= Alpha;
				Total += EdgeToVar_DistEst[i][ EdgeToVar_FindK[i][j] ];
			}
	
			// renormalize for good measure (shouldn't be necessary, could try to take it out)
			for(int j=0;j<EdgeToVar_DistEst[i].size();j++)
			{
				EdgeToVar_DistEst[i][j] *= (1.0/Total);
			}
		}
		else
		{
			EdgeToVar_DistEst[i] = EdgeToVar_DistEst_old[i];
			Alpha = 1.0;
			EdgeToVar_Alpha[i] = Alpha;
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
void KProjectionBP::UpdateBeliefFromVar_Edge(int i)
{
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
			{	
#ifndef NOCHECKS
				{
					vector< FloatType > const & Numbers = EdgeToVar_DistEst[ Var_Fn_Edge[ VarOfEdge[i] ][j] ];
					for(int i=0;i<Numbers.size();i++)
					{
						assert(!isinf(Numbers[i]));
						assert(!isnan(Numbers[i]));
					}
				}
#endif
				EdgeToFn_NonNormal[ i ] [ k ] *=
					EdgeToVar_DistEst[ Var_Fn_Edge[ VarOfEdge[i] ][j] ] [k];
			}
	}

#ifndef NOCHECKS
	{
		vector< FloatType > const & Numbers = EdgeToFn_NonNormal[i];
		for(int i=0;i<Numbers.size();i++)
			assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
	}
#endif
	// normalize the EdgeToFn_Dist
	FloatType MsgTotal = 0.0;
	for(int j=0;j<EdgeToFn_NonNormal[i].size();j++)
	{
		MsgTotal += EdgeToFn_NonNormal[i][j];
	}
	MsgTotal = 1.0 / MsgTotal;
	EdgeToFn_Dist[i] = EdgeToFn_NonNormal[i];
#ifndef NOCHECKS
	{
		vector< FloatType > const & Numbers = EdgeToFn_Dist[i];
		for(int i=0;i<Numbers.size();i++)
			assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
	}
#endif
	for(int j=0;j<EdgeToFn_Dist[i].size();j++)
	{
		EdgeToFn_Dist[i][j] *= MsgTotal;
	}
#ifndef NOCHECKS
	{
		vector< FloatType > const & Numbers = EdgeToFn_Dist[i];
		for(int i=0;i<Numbers.size();i++)
			assert(!isnan(Numbers[i])&&!isinf(Numbers[i]));
	}
#endif
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
vector<FloatType> KProjectionBP::MarginalOfVar(int WhichVar)
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






void KProjectionBP::MarginalOfVars(vector< vector<FloatType> > & Result)
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





void KProjectionBP::DecisionOfVars(vector< int > & Result)
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





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// p to j = (return value) are >= jth value
// j+1 to r are < return value
// int KProjectionBP::partition_largest(vector< FloatType > const & input, vector< int > & Scratch, int p, int r)
// {
// 	FloatType pivot = input[Scratch[p]];

// 	int StoreIndex = r;
// 	int CheckIndex = p;
// 	// everything less or equal to StoreIndex is uncertain
// 	// everything above StoreIndex is less than pivot
// 	while (CheckIndex <= StoreIndex)
// 	{
// 		// if element at CheckIndex less than pivot, put it at StoreIndex
// 		// element now at CheckIndex, was at StoreIndex is still uncertain
// 		// move down StoreIndex
// 		if(input[Scratch[CheckIndex]] < pivot)
// 		{
// 			int tmp = Scratch[CheckIndex];
// 			Scratch[CheckIndex] = Scratch[StoreIndex];
// 			Scratch[StoreIndex] = tmp;
// 			StoreIndex--;
// 		}
// 		// element at CheckIndex can stay put, so move to next
// 		else CheckIndex++;

// 		// eventually CheckIndex == StoreIndex
// 		// if the value there < pivot, then StoreIndex goes down
// 		// otherwise, CheckIndex goes up
// 		// either way, end up at CheckIndex == StoreIndex+1
// 	}

// 	int tmp = Scratch[p];
// 	Scratch[p] = Scratch[StoreIndex];
// 	Scratch[StoreIndex] = tmp;

//   //   for(int i=p;i<=r;i++)
//   //   {
//   //   	if(i<=StoreIndex)
// 		// {
//   //   		if( !(input[Scratch[i]] >= input[Scratch[StoreIndex]]) )
//   //   		{
//   //   			cerr << p << endl;
//   //   			cerr << r << endl;
//   //   			cerr << CheckIndex << endl;
//   //   			cerr << StoreIndex << endl;
//   //   			cerr << pivot << endl;
//   //   			cerr << input[Scratch[StoreIndex]] << endl;
//   //   			cerr << input[Scratch[i]] << endl;
// 	 //    		assert(input[Scratch[i]] >= input[Scratch[StoreIndex]]);
//   //   		}
//   //   	}
//   //   	else
//   //   		assert(input[Scratch[i]] < input[Scratch[StoreIndex]]);
//   //   }

// 	return StoreIndex;
// }

void KProjectionBP::quick_select_largest(vector< FloatType > const & input, vector< int > & Scratch, int p, int r, int k)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if start == end, then there's only one element
	if ( p == r ) return;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// j is the position of the pivot element
	// p to j are greater or equal to the pivot
	// j+1 to r are all strictly less than pivot
	int j = partition_largest(input,Scratch, p, r);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// size of the lesser (left) partition
	// again, j is in the left partition
	int length = j - p + 1;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if the left partition has exactly k elements
	//    then the pivot must be the kth largest value
	// again, the pivot is at j
	if ( length == k ) return;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if k is smaller than the left partition, then the kth largest must be in that partition
	//    furthermore, it is the kth largest from the left partition
	// also, leave out the current pivot.
	else if ( k < length ) quick_select_largest(input,Scratch, p, j - 1, k);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if k is larger than the left partition's size, then the kth largest must be in the right partition
	//    furthermore, leaving off the left block, it is the (k-length)th largest from the right partition
	// also, leave out the current pivot.
	else quick_select_largest(input,Scratch, j + 1, r, k - length);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int KProjectionBP::partition_largest(vector< FloatType > const & input, vector< int > & Scratch, int p, int r)
{
    // FloatType pivot = input[Scratch[r]];
    FloatType pivot = input[Scratch[p]];
    // int p0=p,r0=r;
    
    while ( p < r )
    {
        while ( input[Scratch[p]] > pivot )
            p++;
        
        while ( input[Scratch[r]] < pivot )
            r--;
        
        if ( input[Scratch[p]] == input[Scratch[r]] )
            p++;
        else if ( p < r ) {
            int tmp = Scratch[p];
            Scratch[p] = Scratch[r];
            Scratch[r] = tmp;
        }
    }

    // for(int i=p0;i<=r0;i++)
    // {
    // 	if(i<=r) assert(input[Scratch[i]] >= input[Scratch[r]]);
    // 	else assert(input[Scratch[i]] < input[Scratch[r]]);
    // }
    
    return r;
}

// void KProjectionBP::quick_select_largest(vector< FloatType > const & input, vector< int > & Scratch, int p, int r, int k)
// {
// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// if start == end, then there's only one element
// 	if ( p == r ) return;

// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// j is the position of the pivot element
// 	// p to j are less or equal to the pivot
// 	// j+1 to r are all strictly greater than pivot
// 	int j = partition_largest(input,Scratch, p, r);

// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// size of the lesser (left) partition
// 	// again, j is in the left partition
// 	int length = j - p + 1;

// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// if the left partition has exactly k elements
// 	//    then the pivot must be the kth smallest value
// 	// again, the pivot is at j
// 	if ( length == k ) return;

// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// if k is smaller than the left partition, then the kth smallest must be in that partition
// 	//    furthermore, it is the kth smallest from the left partition
// 	// also, leave out the current pivot.
// 	else if ( k < length ) quick_select_largest(input,Scratch, p, j - 1, k);

// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// if k is larger than the left partition's size, then the kth smallest must be in the right partition
// 	//    furthermore, leaving off the left block, it is the (k-length)th smallest from the right partition
// 	// also, leave out the current pivot.
// 	else quick_select_largest(input,Scratch, j + 1, r, k - length);
// }







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










