//  main.cpp
//  motif prediction
//
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "Node.h"
//#include "Aggregator.h"
#include "myTools.h"
#include<map>
using namespace std;
ifstream infile;
ifstream nodefile;
void readEdge();
void readTime();
int readChosenNodeLine();
int v_current;
int u_current;
int t_current;
int nStates;
int c;
int edgeCount = 1;
int delta = 1000;
int s = 1.5;
int initTime=0;
int timeBlockLength;

int predictStart;
int predictEnd;
int currentReportTime;
bool inPredictionInterval = false;
bool predictionIntervalDone = false;
//tracking multiple overlapping aggregators
std::unordered_map<int,Aggregator*> aggregatorsFinished;
int aggregatorBlockInterval=0;//offset for where the aggregators start
void makePrediction(int t_current, int t);

string predictedOutfileName;
ofstream predictionWindowFile;
ofstream predictionMatrixFile;

int outFileType = -1; //-1 = none, 0 = states, 1 = rates, 2=synthetic
string outFileName; //for outFileType
string motifOutfileName;
std::vector<double> stateSet;
vector<Node*> nodeSet;
int main(int argc, const char * argv[]) {
   if (argc< 8)
   {
	cout<<"ERR, usage: timeLimit edgeCount infile  s  outfile outType OPT:nodeselection"<<endl;
	return 0;
   } 
   int timeLimit = atoi(argv[1]);
   int edgeCount = atoi(argv[2]);
    infile.open(argv[3]);

   c = atof(argv[4]);
   //outFileName = argv[5];
  
   //arg9  = node selected
   std::stringstream ss;
   ss<<argv[5];//<<predictStart<<"to"<<predictEnd<<".txt";
   predictedOutfileName = ss.str();
   std::stringstream ss2;
   ss2<<argv[6];
   motifOutfileName=ss2.str();
   delta=atoi(argv[7]);
   int startmultiplier = atoi(argv[8]);
   int endmultiplier = atoi(argv[9]);
   float blockLengthFrac = atof(argv[10]);
   initTime = atoi(argv[11]);
   timeBlockLength = (int)(blockLengthFrac*(float)timeLimit/2.0)*2;//need to round to even for half offset intervals
   aggregatorBlockInterval = timeBlockLength/2;
   currentReportTime = max(initTime+aggregatorBlockInterval,predictStart-2*aggregatorBlockInterval);
   predictStart = initTime+timeBlockLength*startmultiplier;
   predictEnd  = initTime + timeBlockLength*endmultiplier;
	   
    stateSet = makeStateSet(timeLimit, c, edgeCount);
     for(int i =0; i<stateSet.size();i++){  cout<<"state: "<<i<<" val: "<<stateSet[i]<<endl;}
    nStates = stateSet.size();
    cout<<"Time range: "<<timeLimit<<" starts at "<<initTime<<" timeBlockLength: "<<timeBlockLength<<" edgeCount: "<<edgeCount<<" on file: "<<argv[2]<<" states: "<<nStates<<" delta: "<<delta<<" predicting on " <<predictStart<<" to "<<predictEnd<<" written to "<<predictedOutfileName<<endl;
    std::ofstream myfile;

	Aggregator *A = new Aggregator(stateSet);
	if(outFileType == 2)
	{
		aggregatorBlockInterval = timeBlockLength;
	}
   int edgeIndex =0;
   while( edgeIndex<edgeCount)
   {
		readEdge();

		int maxNew  = max(v_current,u_current);
		//fill in up to what we need- assumes node indices are (almost) sequential from 1
		if(nodeSet.size()<=maxNew)
		{
			for(int i = nodeSet.size();i<=maxNew;i++)
			{
				Node *newv_ptr = new Node(i,stateSet,outFileType,A, timeBlockLength,initTime);
				nodeSet.push_back(newv_ptr);
			}
		}
	  Node *v_vec = nodeSet.at(v_current);
	  Node *u_vec = nodeSet.at(u_current);
	   if(v_vec->getID()!=v_current) cout<<"v err!"<<endl;
	   if (u_vec->getID()!=u_current) cout<<"u err!"<<endl;
	   
	   v_vec->addTimeForBlockCounting(t_current,0,u_current);
	   u_vec->addTimeForBlockCounting(t_current,1,v_current);  
	   while (t_current >= currentReportTime)	   
	   {
		cout<<"reporting rates at time t="<<t_current<<" for an aggregator@ "<<currentReportTime<<endl; 
		 Aggregator *A1 = new Aggregator(stateSet);
		for(int i =0; i<nodeSet.size();i++)
		{
			nodeSet.at(i)->reportRate(predictEnd,timeBlockLength,A1);
		}
		A1->printStateCountCols(t_current);
		aggregatorsFinished.insert(std::make_pair<int&,Aggregator*&>(currentReportTime,A1));	

		cout<<"main calling aggregator to report to synthesizer?"<<outFileType<<","<<t_current<<","<<predictStart<<","<<predictEnd<<endl;


		currentReportTime+=aggregatorBlockInterval;
	  }


	  if(t_current >= predictStart &&!inPredictionInterval)// predictStart>-1)
	  {
		
		//predictStart = -1*predictStart;
		inPredictionInterval = true;
		//if (predictStart ==0) predictStart = -1;
		predictionWindowFile.open(predictedOutfileName,std::ofstream::out);

		cout<<"starting prediction for ("<<predictStart<<","<<predictEnd<<")"<<endl;
		
	  }
	  
	 if(t_current>= predictEnd && !predictionIntervalDone)//end of output interval
	  {

		makePrediction(t_current,edgeIndex);
		edgeIndex = edgeCount+1;//abort loop  
	 }
	edgeIndex++;	  
     }
     if (!predictionIntervalDone)
     {
		cout<<"didn't quite finish prediction interval? "<<edgeIndex<<" "<<t_current<<" for "<<predictStart<<","<<predictEnd<<endl;
		makePrediction(t_current,edgeIndex);
     }



   myfile.close();
    
    
}

void readEdge()
{
    int a = -1;
    int b = -1;
    int t = -1;
    while (a==b && a!=-2)
    {
	a = -2;
        string line;
        getline(infile, line);
        istringstream iss(line);        
        if (!(iss >> a >> b >> t)) { cout<<"NO more edges!\n"; } // error
	
    }
    //cout<<"a :"<<a<<" b: "<<b<<" t: "<<t<<endl;
    v_current = a;
    u_current = b;
    t_current = t;
 	//started and still in interval
    if(inPredictionInterval && t<predictEnd){
	predictionWindowFile<<a<<" "<<b<<" "<<t<<endl;	
    }    
}
int readChosenNodeLine()
{
	 int nid = -1;
	 string line;
	 getline(nodefile,line);
	 istringstream iss(line);
	 if(!(iss>>nid)){ return -1;}
	 

	 return nid;


}
void readTime()
{
    	int t = -1;
	string line;
	getline(infile, line);
	istringstream iss(line);        
	if (!(iss >> t)) { cout<<"NO more edges!\n"; } // error
    	t_current = t;
    
    
}

void makePrediction(int t_current, int t)
{

  


	//predictEnd = -1*predictEnd;
	predictionIntervalDone = true;
	predictionWindowFile.close();
	cout<<"making prediction: "<<nodeSet.size()<<" nodes, time "<<t_current<<endl;


	cout<<"Making (block) prediction at t="<<t<<" for "<<predictStart<<","<<-predictEnd<<endl;
	int lastKey = 0;
	double totalMatrix[6][6];
	for (int i =0; i<6;i++){ for (int j =0; j<6;j++){ totalMatrix[i][j] = 0;}}
	vector<Aggregator *> midAggs;
	vector<Aggregator *> splitAggs;
	int intervalLength = predictEnd -predictStart;
	//int intEnd = -1*predictEnd;
	for(auto kvpair :aggregatorsFinished ) {
		cout<<"Agg key: "<<kvpair.first<<endl;
		int keyTime = kvpair.first;
		if(keyTime<=predictEnd && (keyTime>=(predictStart)+timeBlockLength|| inPredictionInterval))
		{
			cout<<"keytime: "<<keyTime<<"end: "<<predictEnd<<"start: "<<predictStart<<"len: "<<timeBlockLength<<endl;
			if (keyTime==predictEnd && predictStart+timeBlockLength == predictEnd)
			{
				//in split to count 1x
				cout<<kvpair.first<<" only interval; to split"<<endl;
				splitAggs.push_back(kvpair.second);
			}
			else if((predictEnd - keyTime) % timeBlockLength == 0)
			{
				//cout<<kvpair.first<<" to midAggs"<<endl;
				//midAggs.push_back(kvpair.second);		
			}
			else if(kvpair.first > (predictStart)+timeBlockLength)//don't want split @ 1/2 of first interval
			{
				//cout<<kvpair.first<<" to splitAggs"<<endl;
			//	splitAggs.push_back(kvpair.second);
			}
		}
	}
	cout<<"# of each agg type: split: "<<splitAggs.size()<<" mid: "<<midAggs.size()<<endl;
	for(int i =0; i<midAggs.size();i++)
	{
//		double predMatrixLong[6][6];
		double predMatrixReg[6][6];
		Aggregator *midAgg = midAggs.at(i);
		
		//endpoint
		cout<<"makeMotifPrediction: for int= "<<timeBlockLength<<endl;
		midAgg->makePredictions(delta,0,timeBlockLength,predMatrixReg,1);
		cout<<"printing pred from main: (for "<<timeBlockLength<<")"<<endl;
		for(int i =0; i<6;i++)
		{
			for(int j =0; j<6; j++)
			{
				//totalMatrix[i][j]+=predMatrixLong[i][j];
				totalMatrix[i][j]+=predMatrixReg[i][j];
			}
		}

	}
	for(int i =0; i<splitAggs.size();i++)
	{
		Aggregator *splitAgg = splitAggs.at(i);
		double predMatrixReg[6][6];
		double predMatrixHalf[6][6];
		cout<<"makeMotifPrediction: for int(split)="<<timeBlockLength<<endl;
		splitAgg->makePredictions(delta,0,timeBlockLength,predMatrixReg,1);
	
		if(splitAggs.size()>1 || midAggs.size()>0){
			cout<<"makingMotifPrediction: for int/2 (split)"<<timeBlockLength/2<<endl;
			splitAgg->makePredictions(delta,0,timeBlockLength/2,predMatrixHalf,0);
		}
		for(int i =0; i<6;i++)
		{
			for(int j =0; j<6; j++)
			{
				totalMatrix[i][j]+=predMatrixReg[i][j];
				if(splitAggs.size()>1 || midAggs.size()>0){
					totalMatrix[i][j]-=2*predMatrixHalf[i][j];
				}
			}
		}

	}
	predictionMatrixFile.open(motifOutfileName,std::ofstream::out);
	for(int i =0; i<6;i++)
	{
		cout<<"totalpredictrow ";
		for(int j =0; j<6; j++)
		{
			predictionMatrixFile<<totalMatrix[i][j]<<" ";
			cout<<totalMatrix[i][j]<<" ";
		}
		predictionMatrixFile<<endl;
		cout<<endl;
	}
	predictionMatrixFile.close();


}
