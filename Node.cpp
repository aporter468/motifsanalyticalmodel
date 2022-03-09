#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <limits.h>
#include "Node.h"
#include <unordered_map>

using namespace std;
Node::Node() : nodeId(0) {
    outEdgeCount = 0;
    inEdgeCount = 0;
}




Node::Node(int idIn, std::vector<double> states, int outputTypeIn, Aggregator* Ain, int timeBlockLengthIn, int initTime):nodeId(idIn)
{
	A = Ain;
	outEdgeCount = 0;
	inEdgeCount = 0;
	
	timeBlockLength = timeBlockLengthIn;
	lastBlockEnd = initTime+timeBlockLength;
	outputType = outputTypeIn;
        stateSet = states;
	numStates = stateSet.size();
	for (int i =0; i<2; i++)
	{
		lastBlockRates.push_back(-1);
		timeBlockCount[i] = 0;
		
		vector<int> stateset;
		optStates.push_back(stateset);
		vector<double> probset;
		optProbs.push_back(probset);
		vector<int> etimes;
		edgeTimes.push_back(etimes);
		vector<int> egaps;
		edgeGaps.push_back(egaps);
		for (int j = 0; j<numStates; j++)
		{
			optStates[i].push_back(-1);
			optProbs[i].push_back(1.0);
		}

		lastBestState.push_back(-1);
		prevGaps.push_back(-1);
	

		currentOptState[i] = 0;
		currentOptMean[i] = states.at(0);
		sameTime[i] = -1;
		sameTimeCounter[i] = -1;
		firstReportTimePrev[i] = -1;

	}

}


Node::Node(int idIn, std::vector<double> states, int outputTypeIn):nodeId(idIn)
{
   Aggregator *Anew = new Aggregator(states);
    Node(idIn,states,outputTypeIn,Anew,0,100);

}
void Node::addTimeForBlockCounting(int t, int dir, int otherEndID)
{
	if (dir ==0)//is source
	{
		edgeDests.push_back(otherEndID);
	}
	if (timeBlockCount[dir]==0)
	{
		firstReportTime[dir] = t;
	}
	lastReportTime[dir] = t;
	timeBlockCount[dir]++;

}

int Node::getID()
{
	return nodeId;
}
double Node::distributionProbability(double stateP,double lastGap)
{
	double useGap = lastGap;
	double rt = stateP*useGap;
	double prob = exp(-1*rt)*rt;
	return prob;
}

void Node::reportRate(int t, int recordingTime, Aggregator *Ain)
{
		float countSum0 = timeBlockCount[0];
		float countSum1 = timeBlockCount[1];
		
		float elapsed0 = recordingTime;
		float elapsed1 = recordingTime;

		if(oldCounts0.size()>0)	{
			countSum0+=oldCounts0[oldCounts0.size()-1];
		}
		if(oldCounts1.size()>0)	{
			countSum1+=oldCounts1[oldCounts1.size()-1];
		}
	
		float blockRate0 = countSum0/(float)elapsed0;
		float blockRate1 = countSum1/(float)elapsed1;
		vector<float> newRatePair;
		newRatePair.push_back(blockRate0); newRatePair.push_back(blockRate1);
		cout<<"report "<<blockRate0<<","<<blockRate1<<" instead of "<<countSum0/(float)recordingTime<<","<<countSum1/(float)recordingTime<<endl;
		statesFromAgg = Ain->reportOneTimeState(newRatePair,timeBlockCount[0],timeBlockCount[1],t,nodeId);
		
		oldCounts0.push_back(timeBlockCount[0]);
		oldCounts1.push_back(timeBlockCount[1]);

		timeBlockCount[0] = 0;
		timeBlockCount[1] = 0;
		firstReportTimePrev[0] = firstReportTime[0];
		firstReportTimePrev[1] = firstReportTime[1];
		firstReportTime[0] = -1;
		firstReportTime[1] = -1;
}
void Node::reportEdges(Aggregator* Ain)
{
	for (int i =0; i<edgeDests.size();i++)
	{
		Ain->reportEdge(nodeId,edgeDests.at(i));	
	}
	for(int i =0; i<prevEdgeDests.size();i++)
	{
		Ain->reportEdge(nodeId,prevEdgeDests.at(i));
	}
	prevEdgeDests = edgeDests;
	edgeDests.clear();
}

int Node::getCurrentState(int dir)
{
	return lastBestState[dir];
}

int Node::getCurrentStateFromAgg(int dir)
{
	return statesFromAgg[dir];
}
