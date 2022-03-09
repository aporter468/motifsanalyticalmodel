#include <unordered_map>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "Aggregator.h"
class Node
{
    
public:
    Node();

    Node(int idIn, std::vector<double> states, int outputTypeIn);
    Node(int idIn, std::vector<double> states,int outputTypeIn, Aggregator* Ain,  int timeBlockLengthIn, int initTime);
    void addTimeForBlockCounting(int t,int dir, int otherEndID);
    int getID();
    void reportRate(int t, int recordingTime, Aggregator* Ain);
    int getCurrentState(int dir);
    int getCurrentStateFromAgg(int dir);
    void reportEdges(Aggregator* Ain);
private:
   int outEdgeCount;
   int inEdgeCount;
   int nodeId;
   int numStates;
   int currentOptState[2];
   int firstTime = -1;
   int outputType = -1;
   int outputStateSelected = 19;
   int rateWindowSize = 0;
   int sameTimeCounter[2];
   int sameTime[2];
   int timeBlockCount[2];
   int firstReportTime[2];
   int firstReportTimePrev[2];
   int lastReportTime[2];
   std::vector<float> lastBlockRates;
   int timeBlockLength;
   int lastBlockEnd;
   double currentOptMean[2];
   Aggregator *A;
   std::vector<int> lastBestState;
   std::vector<double> stateSet;
   std::vector<std::vector<int> > edgeTimes; 
   std::vector<std::vector<int> > edgeGaps;
   std::vector<std::vector<double> > optProbs;
   std::vector<std::vector<int> > optStates;
   std::vector<int> statesFromAgg;
   std::vector<double > rateWindow;
   std::vector<double> prevGaps;
   std::vector<int> oldCounts0;
   std::vector<int> oldCounts1; 
   std::vector<int> edgeDests;
   std::vector<int> prevEdgeDests;
   double distributionProbability(double stateP,double lastGap);

}; 

