#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <math.h>
using namespace std;

static void getMotifPredictionMatrixBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates, int z, int n, int t,int delta, double (&predictions)[6][6] );
static double computeExpectedMotifsBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates,int k, int z, int n, int t,int delta, int motif[3][2],vector<float> &comboCounts);
static double computeExpectedMotifsSubDeltaBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates,int k, int z, int n, int t,int delta, int motif[3][2], vector<vector<int> > chiset,vector<vector<vector<int> > > sChi,vector<float> &comboCounts);

static void printMotifMatrix(double (&predictions)[6][6]);
static int chooseFunction(int n, int k);


static void makeChiSet(int c, int k, vector<vector<int> > &chiset,vector< vector<vector<int> > > &sChi);
static float getComboCountSet(int c, int chicount, vector<vector<int> > &cSizes, vector<vector<vector<int> > > &sChi,vector<float> &comboCounts);



static vector<double> makeStateSet(int T, int c, int E)
{

    vector<double> states;
   states.push_back(0.0);
  int powChoices[] = {1, 8,4,6,2,5,3,7};

  for(int i =0; i<c; i++)
  {	
	int exponent=i;
	if(i<8)
		exponent =powChoices[i];
	
	float rateval = 1./pow(10,exponent);
	states.push_back(rateval);	
  }

   sort(states.begin(), states.end());

return states;

}


static long int factorial(int n)
{
	if (n<=0)
		return 1;
	if (n<=2)
		return n;
	else
		return n*factorial(n-1);
}
static int chooseFunction(int n, int k)
{
 	if (n<=0 or k<=0)
		return 1;
	if (k>=n)
		return 1;

	int prod = 1;
	for (int i =0; i<k; i++)
	{
		prod*= (n-i);
	}	

	int denom = factorial(k);
	
	int result= (int)((double)prod/(double)denom);
	return result;
}

static void makeChiSet(int c, int k, vector<vector<int> > &chiset,vector< vector<vector<int> > > &sChi)
{
	int k2 = 2*k;
	int chicount = (int)pow(c,k2);
	cout<<"chi count: "<<chicount<<endl;
	for(int i =0; i<chicount;i++)
	{
		int digits[k2];
		int itrack = i;
		vector<int> chival(k2,0);
		std::vector<vector<int> > schi;
		for(int i=0; i<c;i++){
			vector<int> schirow(c,0);
			schi.push_back(schirow);
		}
		for(int d = k2-1; d>-1;d--)//digit
		{
		    digits[d] = (int)(floor(itrack/(pow(c,d))));
		    itrack-=digits[d]*pow(c,d);
		    chival[d] = digits[d];
		  
		}
		for(int d = 0; d<k;d++)
		{
		  schi[chival[d]][chival[d+k]]++;//out groups@1...d, in groups@ d+1...d+k
		}	
		chiset.push_back(chival);
		sChi.push_back(schi);
	}


}


static float getComboCountSet(int c,int chicount, vector<vector<int> > &cSizes, vector<vector<vector<int> > > &sChi,vector<float> &comboCounts)
{
	for(int i =0; i<chicount;i++)
	{
	
		float combocount = 1.0;
			for (int c1 = 0; c1<c; c1++)
			{
				for(int c2 =0; c2<c; c2++)
				{
					int currentChiSize = sChi[i][c1][c2];
					int currentSetSize = cSizes[c1][c2];
					if(currentChiSize > 0 && currentSetSize>0)
					{
						float choose =(float)factorial(currentChiSize)*
						(float)chooseFunction(currentSetSize,currentChiSize);
						combocount= combocount*choose;
					}
					else if (currentSetSize ==0 && currentChiSize>0){
						combocount = combocount*0;
					}
				
				}
			}
		comboCounts[i] = combocount;
	}
}
static double computeExpectedMotifsBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates,int k, int z, int n, int t,int delta, int motif[3][2], vector< vector< int> > chiset, vector< vector<vector< int> > > sChi,vector<float> &comboCounts)
{

	cout<<"compute expected-bidirection: "<<endl;
	cout<<cSizes.size()<<endl;
	int c = cSizes[0].size();
	cout<<"c: "<<c<<"k: "<<k<<endl;
	int k2 = 2*k;
	int chicount = (int)pow(c,k2);


	if(t<=delta){
	 cout<<"t: "<<t<<" delta: "<<delta<<" return subdelt"<<endl;
	 return computeExpectedMotifsSubDeltaBi(cSizes,rates,k,z,n,t,delta,motif,chiset,sChi,comboCounts);
	}
	double subDelta = computeExpectedMotifsSubDeltaBi(cSizes,rates,k,z,n,delta,delta,motif,chiset,sChi,comboCounts);
	cout<<"t>delta subdelta outcome: "<<subDelta<<endl;
	float diff = 0.0;
	int dir = 0;
	cout<<"looping over chivals"<<endl;
	for (int i =0; i<chicount; i++)
	{

		float combocount = comboCounts[i];
		for (int j =0; j<z;j++)//all rates, then by different time windows afterward
		{

			int u = motif[j][0];
			int v = motif[j][1];
			int uC = chiset[i][u];
			int vC = chiset[i][v+k];
			float zrate = rates[uC][vC];//int 0 to t of rate
			combocount = combocount*zrate;
		}
		combocount = combocount*(float)pow(delta,(z-1))*(float)(t-delta);
		diff = diff+combocount;

	}
	diff = diff*0.5;//probability of ordering factor for z =3 

	return subDelta+diff;
}

static double computeExpectedMotifsSubDeltaBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates,int k, int z, int n, int t,int delta, int motif[3][2], vector<vector<int> > chiset,vector<vector<vector<int> > > sChi, vector<float> &comboCounts)
{
	printf("Called Compute Expected subdelta (Bi): ");
	for (int i =0; i<3;i++)
	{
		printf("(%d,%d) ",motif[i][0],motif[i][1]);
	}
	printf(" k  = %d\n",k);
	double expectedCount = 1;
	int c = cSizes[0].size();
	int k2 = 2*k;
	printf("c: %d k: %d k2: %d\n",c,k,k2);
	int chicount = (int)pow(c,k2);
	printf("# chi values (c^2k): %d chisetlen %d\n",chicount,chiset.size());
	double countsum = 0.0;
	for (int i =0; i<chicount;i++)
	{
		double motifprob = 1.0;

		if(comboCounts[i]>0)
		{
			for (int j = 0; j<z;j++)
			{
				int u= 0;
				int v = 0;
				u = motif[j][0];
				v = motif[j][1];

				int uC = chiset[i][u];
				int vC = chiset[i][k+v];//second half of chi is "in"
				float zrate = rates[uC][vC];//int 0 to t of rate
				motifprob *=(zrate*t);
			}
			float orderprob = 1.0/6.0;//always the same for constant rate functions
			float combocount = motifprob*orderprob*comboCounts[i];
			countsum+=combocount;
		}	
	}
	printf("predicted count: %f\n",countsum);


	return countsum;


}


 

static void getMotifPredictionMatrixBi(vector<vector<int> > &cSizes, vector<vector<double> > &rates, int z, int n, int t,int delta, double (&predictions)[6][6] )
{
	cout<<"getMotifPredictionMatrixBi"<<endl;
	int edge1[2] ={0,1};
	int edgeList[6][2] = {{2,1},{1,2},{2,0},{0,2},{1,0},{0,1}};
	//double predictions[6][6];
	vector<vector<int> > chiset2;
	vector<vector<vector<int> > > sChi2;
	int c = cSizes[0].size();
	int totalN = 0;
	cout<<"t: "<<t<<" delta: "<<delta<<endl;
	cout<<"rates recvd: "<<endl;
	for(int i =0; i<c; i++)
	{
		for(int j =0; j<c;j++)
		{
			cout<<rates[i][j]<<" ";

		}
		cout<<endl;
	}
	cout<<"Category sizes:"<<endl;
	for(int i =0; i<c;i++)
	{
		for(int j =0; j<c; j++)
		{
		cout<<cSizes[i][j]<<" ";
		totalN+=cSizes[i][j];
		}
		cout<<endl;
	}
	cout<<"total nodes counted: "<<totalN<<endl;


	makeChiSet(c,2,chiset2,sChi2);
	vector<vector<int> > chiset3;
	vector<vector<vector<int> > >sChi3;
	makeChiSet(c,3,chiset3,sChi3);
 	vector<float> comboCounts2(chiset2.size(),0);
	vector<float> comboCounts3(chiset3.size(),0);
	getComboCountSet(c,chiset2.size(), cSizes,sChi2,comboCounts2);
	cout<<"got 2 combos"<<endl;
	getComboCountSet(c,chiset3.size(),cSizes,sChi3,comboCounts3);
	cout<<"got 3 combos"<<endl;
	for(int twoIter = 0; twoIter<6; twoIter++)
	{
		for (int threeIter = 5; threeIter>=0; threeIter--)
		 {
			int motif[3][2];
			motif[0][0] = edge1[0]; motif[0][1] = edge1[1];
			motif[1][0] = edgeList[twoIter][0]; motif[1][1] = edgeList[twoIter][1];
			motif[2][0] = edgeList[threeIter][0]; motif[2][1] = edgeList[threeIter][1];
			int k =2;
			if (motif[1][0]==2 || motif[1][1] ==2 || motif[2][0]==2 || motif[2][1]==2)
			{
				k = 3;	
			}
			cout<<"motif: "<<twoIter<<","<<(5-threeIter)<<endl;
			double predict = 0;
			if (k==2){
				 predict = computeExpectedMotifsBi(cSizes,rates,k,z,n,t,delta,motif,chiset2,sChi2,comboCounts2); 
			}
			else{
				predict = computeExpectedMotifsBi(cSizes,rates,k,z,n,t,delta,motif,chiset3,sChi3,comboCounts3);
			}
			cout<<"returned predict: "<<predict<<endl;
			 predictions[twoIter][5-threeIter] = predict;	
		}	
	}
	cout<<"finished matrix: "<<endl;
	for(int twoIter = 0; twoIter<6; twoIter++)
	{
		cout<<"motifrow ";
		for (int threeIter = 0; threeIter<6; threeIter++)
		 {
			cout<<predictions[twoIter][threeIter]<<" ";
		}
		cout<<endl;
	}	
}
static void printMotifMatrix(double (&predictions)[6][6])
{
	printf("Motif Predictions:\n");
	for (int i = 0; i<6; i++)
	{
		for (int j =0; j<6; j++)
		{
			printf("%f ",predictions[i][j]);
		}
		printf("\n");
	}

}

