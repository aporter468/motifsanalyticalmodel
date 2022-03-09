#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <limits.h>
#include "Aggregator.h"
#include "myTools.h"
#include <unordered_map>


using namespace std;
Aggregator::Aggregator() {
}

Aggregator::Aggregator(std::vector<double> states)
{
        stateSet = states;
	resetStateCounts();	
	resetEdgeCounts();
	cout<<"making col state/costs"<<endl;
	for(int i =0; i<2;i++){
		vector<double> stateCostCol(stateSet.size(),0.0);
		stateDistCols.push_back(stateCostCol);
		vector<int> stateCountCol(stateSet.size(),0);
		stateCountCols.push_back(stateCountCol);
	}
	edgeOutCount = 0;
	edgeInCount = 0;
	numNodes = 0;
	stateChangesByDir = 0;
	currentDistanceTotal = 0.;
	distRunningSum =0.;
        cout<<"init states:" <<endl;
 	for(int i =0; i<stateSet.size();i++)
	{
		for(int j =0;j<stateSet.size();j++){
			cout<<stateCounts[i][j]<<" ";
		}
		cout<<endl;
	}


}
void Aggregator::resetStateCounts()
{
	vector<vector<int> > newStateCounts;
	for (int j=0; j<stateSet.size(); j++)
	{
		vector<int> countDirVec;
		for (int i =0; i<stateSet.size();i++)
		{
			countDirVec.push_back(0);
		}
		newStateCounts.push_back(countDirVec);
	}
	stateCounts = newStateCounts;
}
void Aggregator::resetEdgeCounts()
{
	vector<vector<int> > newEdgeMatrix;
	for (int j=0; j<stateSet.size(); j++)
	{
		vector<int> edgeCountVec;
		for (int i =0; i<stateSet.size();i++)
		{
			edgeCountVec.push_back(0);
		}
		newEdgeMatrix.push_back(edgeCountVec);
	}
	edgeMatrix = newEdgeMatrix;

}

void Aggregator::makePredictions(int delta, int start, int end, double (&predictions)[6][6],int doMatrixReduce)
{
	cout<<"states at prediction time: "<<endl;
	for(int i =0; i<stateSet.size();i++)
	{
		for(int j =0;j<stateSet.size();j++){
			cout<<stateCounts[i][j]<<" ";
		}
		cout<<endl;
	}
	float interval = (float)(end-start);
	if(doMatrixReduce>0)
		reduceGroupMatrix();
	rateMatrix = makeRateMatrix(interval,1);
	printRateMatrix();
	printEdgeMatrix();
	cout<<"Agg make prediction for: "<<start<<","<<end<<": "<<(end-start)<<endl;	
	getMotifPredictionMatrixBi(stateCounts, rateMatrix,3,numNodes,end-start,delta,predictions);
  	printMotifMatrix(predictions);	
	

}
void Aggregator::clearStateCounts()
{
	cout<<"Aggregator clears state counts"<<endl;
	for(int i =0; i<stateSet.size();i++)
	{
		stateCountCols[0][i]=0;
		stateCountCols[1][i]=0;
		stateDistCols[0][i]=0;
		stateDistCols[1][i]=0;
		for(int j =0 ;j<stateSet.size();j++)		
		{
			stateCounts[i][j]=0;
		}
	}	
}

//counting type = 0 -> return rates by observed edges; 1-> return rates by averaging out rates
vector<vector<double> > Aggregator::makeRateMatrix(float interval,int countingType)
{
	vector<vector<double> >stateAvgCosts;
	float dirSums[2];
	float avgSums[2];
	int vertexCounts[2];
	for(int i =0;i<2;i++){
		vector<double> avgsCol;
		dirSums[i] = 0.0;
		avgSums[i]=0.0;
		vertexCounts[i]=0;
		for(int j = 0;j<stateDistCols[i].size();j++)
		{
			float nextAvg = 0.0;
			if(stateCountCols[i][j]>0)
				nextAvg= stateDistCols[i][j] /(float)stateCountCols[i][j];
			cout<<"col: "<<i<<" row: "<<j<<" sum: "<<stateDistCols[i][j]<<" count: "<<stateCountCols[i][j]<<" avg: "<<nextAvg<<endl;
			avgsCol.push_back(nextAvg);
			
			dirSums[i]+=stateDistCols[i][j];
			avgSums[i]+= stateDistCols[i][j]/(float)stateCountCols[i][j];
			vertexCounts[i]+=stateCountCols[i][j];
		}
		stateAvgCosts.push_back(avgsCol);
	}
	cout<<"dir sums: "<<dirSums[0]<<","<<dirSums[1]<<" predictededges: "<<dirSums[0]*interval<<","<<dirSums[1]*interval<<" for interval "<<interval<<" vertex counts: "<<vertexCounts[0]<<","<<vertexCounts[1]<<endl;
	cout<<"pred.edges "<<dirSums[0]*interval<<endl;

	double stateSum = 0.0;
	for(int i =0; i<stateSet.size();i++){
		stateSum+=stateSet[i];
	}


	vector<int> outEdgeSums(stateSet.size(),0);
	vector<int> inEdgeSums(stateSet.size(),0);
	for(int i =0; i<stateSet.size();i++)
	{
		for(int j =0; j<stateSet.size();j++)
		{
			outEdgeSums[i]+=edgeMatrix[i][j];
			inEdgeSums[j] += edgeMatrix[i][j];
		}
	}
	cout<<"rates from states for comparison: "<<endl;
	vector<vector<double> > rateMatrix;
	int edgeMatrixTotal=0;
	for(int i =0; i<stateSet.size();i++)	
	{
		vector<double> rateVec;
		for(int j = 0; j<stateSet.size();j++)
		{
			float totalprod = (float)stateCountCols[0][i]*(float)stateCountCols[1][j];
	
			cout<<"matrix"<<i<<","<<j<<": "<<stateAvgCosts[0][i]<<" "<<stateDistCols[1][j]<<" "<<dirSums[1]<<endl;
			float rateij = stateAvgCosts[0][i]*stateDistCols[1][j]/(dirSums[1]);
			if(stateCountCols[1][j]>0)
			{
				if(i==j)
				{
					rateij = rateij/(stateCountCols[1][j]-1);
					if (stateCountCols[1][j]==1) rateij=0;
				}
				else
				{
					rateij = rateij/stateCountCols[1][j];
				}
			}
			else
			{
				rateij = 0;
			}
			rateVec.push_back(rateij);
			
			edgeMatrixTotal+=edgeMatrix[i][j];
		}
		rateMatrix.push_back(rateVec);
	}

       	cout<<"edge matrix total: "<<edgeMatrixTotal<<endl;
	for(int i = 0; i<stateSet.size();i++)
	{
		cout<<"ratesfromedges: ";
		vector<double> edgeRateRow;
		for(int j = 0; j<stateSet.size();j++)
		{
			float edgeRateij = 0.;
			if(edgeMatrix[i][j]!=0){
				edgeRateij = edgeMatrix[i][j]/(interval*stateCountCols[0][i]*stateCountCols[1][j]);	
			}
			edgeRateRow.push_back(edgeRateij);
			cout<<edgeRateij<<" ";
		}
		edgeRateMatrix.push_back(edgeRateRow);
		cout<<endl;
	}

	float edgePredictionTotal = 0.;
	for(int i =0; i<stateSet.size();i++){
		float colSum = 0.0;
		for(int j =0; j<stateSet.size();j++){
		 colSum+=edgeRateMatrix[i][j]*stateCountCols[0][i]*stateCountCols[1][j];
		}
		float edgepred = interval*colSum;
	
		edgePredictionTotal+=edgepred;
	}
	
	cout<<"predicted total: "<<edgePredictionTotal<<endl;
	if (countingType ==0)
		return edgeRateMatrix;
	else
		return rateMatrix;
}

void Aggregator::changeBlockRateStates(std::vector<float> prevRates, std::vector<float> newRates, int nid, int blockEndTime){
	for(int i = 0; i<2;i++) 
	{
		cout<<i<<"i:"<<nid<<" "<<blockEndTime<<" "<<newRates[i]<<endl;
		if(prevRates[i] ==-1) prevRates[i] = 0;
	}
	
	int prevState[2];
	int newState[2];
	for(int i =0; i<2;i++){
		prevState[i] = getStateFromRate(prevRates[i]);
		newState[i] = getStateFromRate(newRates[i]);
	}
	stateCounts[prevState[0]][prevState[1]]--;
	stateCounts[newState[0]][newState[1]]++;
	


	for(int i =0;i<2;i++){
		stateDistCols[i][prevState[i]]-=prevRates[i];//even if prev not exist ok- would be 0
		stateDistCols[i][newState[i]]+=newRates[i];
		stateCountCols[i][prevState[i]]--;	
		stateCountCols[i][newState[i]]++;
	}
		
}
int Aggregator::getStateFromRate(float rate)
{
  float minDiff = abs(rate - stateSet[0]);
  int bestState =0;
  for(int i =0; i<stateSet.size();i++){
     float newDiff = abs(rate - stateSet[i]);
     if(newDiff<minDiff){
	minDiff = newDiff;
	bestState = i;
     }
  }
  return bestState;
}
void Aggregator::setInitZeroState()
{
	stateCounts[0][0]++;
	stateCountCols[0][0]++;
	stateCountCols[1][0]++;
}
vector<int> Aggregator::reportOneTimeState(vector<float> newRates, int outCount, int inCount, int t, int id)
{
	vector<int> newState;
	
	for(int i =0; i<2;i++){
		newState.push_back(getStateFromRate(newRates[i]));
	}

	vertexStates.insert(make_pair<int&,vector<int>& >(id,newState));
	cout<<"rate reported: "<<newRates[0]<<","<<newRates[1]<<" to states: "<<newState[0]<<","<<newState[1]<<endl;
	stateCountCols[0][newState[0]]++;
	stateCountCols[1][newState[1]]++;
	stateCounts[newState[0]][newState[1]]++;
	stateDistCols[0][newState[0]]+=newRates[0];
	stateDistCols[1][newState[1]]+=newRates[1];
    	edgeOutCount += outCount;
 	edgeInCount+=inCount;
	return newState;

}

void Aggregator::countEdgeEnds(int ustate,int vstate)
{
	cout<<"count edge: s1:"<<ustate<<" s2:"<<vstate<<endl;
	if(ustate <0) ustate = 0;
	if(vstate <0) vstate = 0;
	if(ustate>-1 && vstate>-1)	
		edgeMatrix[ustate][vstate]++;
}

void Aggregator::printEdgeMatrix()
{
	cout<<"Edge Count Matrix: "<<endl;
	for(int i =0; i<stateSet.size();i++)
	{
		for(int j =0; j<stateSet.size();j++)
		{
			cout<<edgeMatrix[i][j]<<" ";
		}
		cout<<endl;
	}

}
void Aggregator::printRateMatrix()
{
	cout<<"Rate Matrix: "<<endl;
	for(int i =0; i<stateSet.size();i++)
	{
		for(int j =0; j<stateSet.size();j++)
		{
			cout<<rateMatrix[i][j]<<" ";
		}
		cout<<endl;
	}

}


void Aggregator::reportEdge(int source, int dest)
{
	int srcState =vertexStates[source][0];
	int destState = vertexStates[dest][1]; 
	countEdgeEnds(srcState, destState);
}

void Aggregator::printStateCountCols(int time)
{
     for(int j =0; j<2;j++){
	for(int i=0; i<stateSet.size();i++)
	{
		float avgStateRate = stateDistCols[j][i]/stateCountCols[j][i];
		cout<<"state"<<i<<"dir"<<j<<"count: "<<time<<" "<<stateCountCols[j][i]<<" "<<avgStateRate<<endl;
	}
     }
}

void Aggregator::reduceGroupMatrix()
{
	cout<<"matrix reductions: "<<endl;
	vector<int> removeGroups;
	vector<int> keepGroups;
	for(int i =0; i<stateCounts.size();i++)
	{
		if (stateCountCols[0][i] == 0 && stateCountCols[1][i] ==0)
		{
			removeGroups.push_back(i);
			cout<<"remove group: "<<i<<endl;
		}
		else
		{
			keepGroups.push_back(i);
		}
	}	
	vector<vector<int> > newGroupMatrix;
	vector<vector< int> >newCountCols;
	vector<vector< double> >newDistCols;
	for(int i =0; i<2;i++)
	{
		vector<int> countCol;
		newCountCols.push_back(countCol);
		vector<double> distCol;
		newDistCols.push_back(distCol);
	}
	for(int i =0; i<keepGroups.size();i++)
	{
		vector<int> newMatrixRow;
		for(int j =0;j <keepGroups.size();j++)	
		{
			newMatrixRow.push_back(stateCounts[keepGroups[i]][keepGroups[j]]);
		}
		newGroupMatrix.push_back(newMatrixRow);
		for(int j =0;j <2;j++)
		{
			newCountCols[j].push_back(stateCountCols[j][keepGroups[i]]);
			newDistCols[j].push_back(stateDistCols[j][keepGroups[i]]);
		}
	}


	cout<<"new matrix:" <<endl;
	for(int i =0; i<keepGroups.size();i++)
	{
		for(int j =0;j <keepGroups.size();j++)	
		{
			cout<<newGroupMatrix[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"new cols: "<<endl;
	for(int i =0;i<keepGroups.size();i++)
	{
		cout<<"cols: "<<newCountCols[0][i]<<","<<newDistCols[0][i]<<"   "<<newCountCols[1][i]<<","<<newDistCols[1][i]<<endl;
	}

	stateCounts = newGroupMatrix;
	stateDistCols = newDistCols;
	stateCountCols = newCountCols;
}
