#include "StochasticBP2.h"

#ifndef NULL
#define NULL 0
#endif

// DONE
// Verified
StochasticBP2::StochasticBP2(FactorGraph newGraph)
{
	UseGraph(newGraph);
}

string StochasticBP2::PrintParams()
{
	return string("none");
}


// DONE
int StochasticBP2::UseGraph(FactorGraph newGraph)
{
	GraphSolver::UseGraph(newGraph);

	int EdgeCount = VarOfEdge.size();

	EdgeToVarMsg.resize(EdgeCount);
	EdgeToFnMsg.resize(EdgeCount);
	EdgeToVarSymbols.resize(EdgeCount);

	int WhichEdge = 0;
	// go through each function node
	for(int i=0;i<graph.functions.size();i++)
	{
		// for every variable connected to that function node
		// i.e. every edge from that function node
		for(int j=0;j<graph.functions[i].vars.size();j++)
		{
			int tempSize;

			// divide the number of fn values by that variable's cardinality
			tempSize = graph.functions[i].ValuesSize /
				graph.variables[ graph.functions[i].vars[j] ].cardinality;

			tempSize = graph.variables[ graph.functions[i].vars[j] ].cardinality;
			// initial assignment of messages is UNIFORM
			EdgeToFnMsg[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) );
			EdgeToVarMsg[WhichEdge].assign( tempSize , 1.0L / (FloatType)(tempSize) );

			EdgeToVarSymbols[WhichEdge].resize(graph.functions[i].vars.size());
			
			WhichEdge++;
		}
	}

	IterNumber = 0;

	return 0;
}




void StochasticBP2::ResetState()
{
	for(auto & Msg:EdgeToVarMsg)
	{
		Msg.assign( Msg.size() , 1.0L / (FloatType)(Msg.size()) );
	}
	IterNumber = 0;
}


void StochasticBP2::GetState(vector< vector< FloatType > > & State)
{
	// State consists of:
	//  - EdgeToVarMsg[][] 
	State.resize(EdgeToVarMsg.size()+1);
	int Index=0;
	for(auto & v:EdgeToVarMsg)
	{
		State[Index++] = v;
	}
	State[Index].push_back(FloatType(IterNumber));
}


void StochasticBP2::RestoreState(vector< vector< FloatType > > const & State)
{
	int Index=0;
	for(auto & v:EdgeToVarMsg)
	{
		v = State[Index++];
	}
	IterNumber = int(State[Index][0]);
}





void StochasticBP2::Iter()
{
	IterNumber++;
	Iter(2.0 / (1.0 + IterNumber));
}


void StochasticBP2::Iter(FloatType StepSize)
{
	// update every EdgeToFnMsg
	// structure variables modified in loop iteration 'i':
	//   - EdgeToFnMsg[i]
	// no other EdgeToFnMsg is accessed in iter 'i'
	// -> loop iterations are independent, can safely multithread
	for(int i=0;i<EdgeToFnMsg.size();i++)
	{
		UpdateMsgFromVar_Edge(i);
	}

	// generate EdgeToVarSymbols by sampling from joint distribution
	// then update the EdgeToVar messages using the step update
	uniform_real_distribution<FloatType> d(0.0,1.0);
	for(int i=0;i<EdgeToVarSymbols.size();i++)
	{
		UpdateMsgFromFn_Edge(i,StepSize,r,&d);
	}
}































///////////////////////
struct ThreadData
{
	StochasticBP2 * Graph;
	int * Progress;
	FloatType StepSize;
	pthread_mutex_t * Mutex;
	default_random_engine * r;
	uniform_real_distribution<FloatType> * d;
};
///////////////////////
int const ChunkSize = 15;
///////////////////////

void * EdgeToFnThread_StochasticBP2(void * d)
{
	// cast the input pointer so we can use it
	ThreadData * theData = (ThreadData *)d;
	
	int QuitThread;
	do
	{
		int StartCompute;
		QuitThread = 0;
		pthread_mutex_lock(theData->Mutex);
		if(*(theData->Progress) < theData->Graph->EdgeToFnMsg.size())
		{
			StartCompute = *(theData->Progress);
			*(theData->Progress) += ChunkSize;
		}
		else QuitThread = 1;
		pthread_mutex_unlock(theData->Mutex);
		if(!QuitThread)
		{
			for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToFnMsg.size();i++)
				theData->Graph->UpdateMsgFromVar_Edge(i);
		}
	}while(!QuitThread);
	return NULL;
}


void * EdgeToVarThread_StochasticBP2(void * d)
{
	// cast the input pointer so we can use it
	ThreadData * theData = (ThreadData *)d;
	
	int QuitThread;
	do
	{
		int StartCompute;
		QuitThread = 0;
		pthread_mutex_lock(theData->Mutex);
		if(*(theData->Progress) < theData->Graph->EdgeToVarMsg.size())
		{
			StartCompute = *(theData->Progress);
			*(theData->Progress) += ChunkSize;
		}
		else QuitThread = 1;
		pthread_mutex_unlock(theData->Mutex);
		if(!QuitThread)
		{
			for(int i=StartCompute;i<StartCompute+ChunkSize && i<theData->Graph->EdgeToVarMsg.size();i++)
				theData->Graph->UpdateMsgFromFn_Edge(i,theData->StepSize,theData->r,theData->d);
		}
	}while(!QuitThread);
	return NULL;
}



void StochasticBP2::IterThreads(int NumThreads)
{
	IterNumber++;
	IterThreads(2.0 / (1.0 + IterNumber) , NumThreads);
}

void StochasticBP2::IterThreads(FloatType StepSize,int NumThreads)
{
	vector< default_random_engine > ThreadRNG;
	ThreadRNG.resize(NumThreads);
	vector< uniform_real_distribution< FloatType > > ThreadDist(NumThreads,uniform_real_distribution< FloatType >(0.0,1.0));
	uniform_int_distribution<long> SeedDist;
	for(int i=0;i<ThreadRNG.size();i++)
	{
		long seed = SeedDist(*r);
		ThreadRNG[i].seed(seed);
	}

	//////////////////////////////////////////////
	// update every EdgeToFnMsg
	{
		vector<pthread_t> myThreads;
		myThreads.resize(NumThreads);

		vector<ThreadData> theData;
		theData.resize(NumThreads);
		pthread_mutex_t Mutex;
		pthread_mutex_init(&Mutex,NULL);
		int Progress=0;

		for(int i=0;i<NumThreads;i++)
		{
			theData[i].Graph = this;
			theData[i].Progress = &Progress;
			theData[i].StepSize = StepSize;
			theData[i].Mutex = &Mutex;
			theData[i].r = &ThreadRNG[i];
			theData[i].d = &ThreadDist[i];
		}

		for(int i=0;i<NumThreads;i++)
		{
			pthread_create( &myThreads[i], NULL, EdgeToFnThread_StochasticBP2, (void*)(&theData[i]) );
		}
		for(int i=0;i<NumThreads;i++)
		{
			pthread_join(myThreads[i],NULL);
		}
		pthread_mutex_destroy(&Mutex);
	}
	//////////////////////////////////////////////
	// update every EdgeToVarMsg
	{	
		vector<pthread_t> myThreads;
		myThreads.resize(NumThreads);

		vector<ThreadData> theData;
		theData.resize(NumThreads);
		pthread_mutex_t Mutex;
		pthread_mutex_init(&Mutex,NULL);
		int Progress=0;

		for(int i=0;i<NumThreads;i++)
		{
			theData[i].Graph = this;
			theData[i].Progress = &Progress;
			theData[i].StepSize = StepSize;
			theData[i].Mutex = &Mutex;
			theData[i].r = &ThreadRNG[i];
			theData[i].d = &ThreadDist[i];
		}

		for(int i=0;i<NumThreads;i++)
		{
			pthread_create( &myThreads[i], NULL, EdgeToVarThread_StochasticBP2, (void*)(&theData[i]) );
		}
		for(int i=0;i<NumThreads;i++)
		{
			pthread_join(myThreads[i],NULL);
		}
		pthread_mutex_destroy(&Mutex);
	}
	//////////////////////////////////////////////
}




























void StochasticBP2::UpdateMsgFromVar_Edge(int i)
{
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




void StochasticBP2::UpdateMsgFromFn_Edge(int i,FloatType StepSize,default_random_engine * r,uniform_real_distribution<FloatType> * d)
{
	// which function does the edge come from?
	int theFn=FnOfEdge[i];
	int LeaveOut = -10000;

	// generate EdgeToVarSymbols by sampling from joint distribution
	// sample independently from the vars attached to the Fn, excluding current edge
	for(int j=0;j<graph.functions[theFn].vars.size();j++)
	{
		// exclude the current edge
		if(Fn_Var_Edge[ theFn ][j] != i)
		{
			FloatType tempSum = 0;
			for(int k=0;k<EdgeToFnMsg[ Fn_Var_Edge[ theFn ][j] ].size();k++)
			{
				tempSum += EdgeToFnMsg[ Fn_Var_Edge[ theFn ][j] ][k];
			}
			FloatType RandU = ((*d)(*r)) * tempSum;
			tempSum = EdgeToFnMsg[ Fn_Var_Edge[ theFn ][j] ][0];
			EdgeToVarSymbols[i][j] = 0;
			while(tempSum < RandU)
			{
				EdgeToVarSymbols[i][j]++;
				tempSum += EdgeToFnMsg[ Fn_Var_Edge[ theFn ][j] ][EdgeToVarSymbols[i][j]];
			}
		}else LeaveOut=j;
	}

	// update the EdgeToVarMsg using the step update
	FloatType MsgTotal = 0.0;
	for(int j=0;j<EdgeToVarMsg[i].size();j++)
	{
		EdgeToVarMsg[i][j] *= (1.0 - StepSize);
		EdgeToVarSymbols[i][LeaveOut]=j;
		EdgeToVarMsg[i][j] += StepSize * graph.values[graph.functions[theFn].values][CoordsToFullIndex(theFn,EdgeToVarSymbols[i])];
		MsgTotal += EdgeToVarMsg[i][j];
	}
	MsgTotal = 1.0 / MsgTotal;
	for(int j=0;j<EdgeToVarMsg[i].size();j++)
	{
		EdgeToVarMsg[i][j] *= MsgTotal;
	}
}
































// DONE
vector<FloatType> StochasticBP2::MarginalOfVar(int WhichVar)
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





void StochasticBP2::MarginalOfVars(vector< vector<FloatType> > & Result)
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




void StochasticBP2::DecisionOfVars(vector< int > & Result)
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










