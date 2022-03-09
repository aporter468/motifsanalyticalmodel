#include <unordered_map>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

class Aggregator
{
    
public:
    Aggregator();
    Aggregator(std::vector<double> states);
    /* index dirs as 0=out, 1 = in */ 
    void changeState(int start, int end, int nid, int dir);
    void makePredictions(int delta, int start,int end, double (&predictions)[6][6],int doMatrixReduce);
    void changeBlockRateStates(std::vector<float> prevRates, std::vector<float> newRates, int nid, int blockEndTime);
    void setInitZeroState();
   void resetStateCounts();
   std::vector<int> reportOneTimeState(std::vector<float> newRates, int outCount, int inCount,int t, int id);   
   void countEdgeEnds(int ustate, int vstate);
   void printEdgeMatrix();

   void printRateMatrix();
   void printStateCountCols(int time);
   void clearStateCounts();
   void reportEdge(int source, int dest);
  

private:
   int numNodes;
   int stateChangesByDir;
   int edgeOutCount;
   int edgeInCount;
   float distRunningSum;
   std::vector<double> stateSet;
   std::vector<std::vector<int> > stateCounts;
   std::vector<std::vector<double> > stateDistCols;
   std::vector<std::vector<int> > stateCountCols;
   std::vector<std::vector<double> > rateMatrix;

   std::vector<std::vector<double> > edgeRateMatrix;
   float currentDistanceTotal;
   std::vector<std::vector<double> >  makeRateMatrix(float interval,int countingType);
   std::vector<std::vector<int> > edgeMatrix;
   int getStateFromRate(float rate);
   std::unordered_map<int, std::vector<int> > vertexStates;
   void resetEdgeCounts();
   void reduceGroupMatrix();
}; 

