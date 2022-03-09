//  main_gen.cpp
//  synthetic temporal network generator
//
//
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "../myTools.h"
#include<map>
using namespace std;

double getNextTime(double lambda);


ifstream rateMatrixFile;
ifstream groupMatrixFile;
ofstream outfile;
vector<vector<double > >rateMatrix;
vector<vector<int> > countMatrix;
vector<int> nodeToInBlock;
vector<int> nodeToOutBlock;
void readMatrices(int groupCount);
int totalNodes;
int T;
int delta;
void sampleGraph();
int main(int argc, const char * argv[]) {
	rateMatrixFile.open(argv[1]);
	groupMatrixFile.open(argv[2]);
	int groupCount = atoi(argv[3]);
	delta = atoi(argv[4]);
	T = atoi(argv[5]);
	outfile.open(argv[6]);
	int testType = atoi(argv[7]);
	readMatrices(groupCount);
	int N = 0;
	srand(time(0));
	for(int i =0; i<groupCount;i++)
	{
		for(int j =0;j<groupCount;j++)
		{
			cout<<rateMatrix[i][j]<<","<<countMatrix[i][j]<<" ";
			N+=countMatrix[i][j];
		}
		cout<<endl;
	}
	if(testType == 0)
	{
		double predictions[6][6];
		cout<<"calling with: N="<<N<<" T="<<T<<" delta="<<delta<<endl;
			getMotifPredictionMatrixBi(countMatrix,rateMatrix,3,N,T,delta,predictions);
		for (int i =0; i<6;i ++)
		{
			for(int j = 0; j<6; j++)	
			{
				outfile<<predictions[i][j]<<" ";
			}
			outfile<<endl;
		}
	}
	else
	{
		sampleGraph();
	}
	outfile.close();
}


void readMatrices(int groupCount)
{
	totalNodes = 0;
	for (int j = 0; j<groupCount; j++)
	{
		vector<int> countRow;	
		vector<double> rateRow;
	
		string line;
		getline(rateMatrixFile,line);
		istringstream iss(line);
		double r;
		while (iss>>r)
		{
			rateRow.push_back(r*10);
		}
		string line2;
		getline(groupMatrixFile,line2);
		istringstream iss2(line2);
		int c;
		while(iss2>>c)
		{
			countRow.push_back(c);
			for (int k = 0; k<c;k++)
			{
				nodeToInBlock.push_back(countRow.size()-1);
				nodeToOutBlock.push_back(countMatrix.size());
			}
			totalNodes +=c;
		}

		rateMatrix.push_back(rateRow);
		countMatrix.push_back(countRow);
	
	}
/*	for (int i =0; i<nodeToOutBlock.size();i++)
	{
		cout<<"out: "<<nodeToOutBlock[i]<<" in: "<<nodeToInBlock[i]<<endl;
	}
*/
}


void sampleGraph()
{
	int subdeltacount = 0;
	for (int i =0; i<totalNodes; i++)
	{
		int iblock = nodeToOutBlock[i];
		for (int j =0; j<totalNodes;j++)
		{
			int jblock = nodeToInBlock[j];
			if (j !=i)
			{

					double P = rateMatrix[iblock][jblock];
					double t = getNextTime(P);
					//printf("%d,%d times: ",i,j);
					int count = 0;

					while (t<T)
					{
					//	printf("%f ",t);
						//temporal_data_[i](j).Add(t);
						outfile<<i<<" "<<j<<" "<<(int)t<<endl;
						double tint =  getNextTime(P);
						if (t<delta){
							subdeltacount++;
						}
						t = t+tint;	
						//count++;
				
					}
					//countfile<<count<<endl;
					//printf("c=%d\n",count);
			}
		}
	}
	printf("edges < delta: %d",subdeltacount);
	outfile.close();

 }
double getNextTime(double lambda)
{
	double u =  (double)rand() / RAND_MAX;	
	double r = -1*log(u)/lambda;
	return r;
}
