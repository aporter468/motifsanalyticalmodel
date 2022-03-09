#generates a synthetic graph given the groups+rates for stochastic block model:
#row corresponds to "from" group and column to "to" group in rates file 
#row corresponds to "out" group and column to "in" group in groups file
import numpy as np
import sys
import random

args = sys.argv
outfileName = str(args[1])
outfile = open(outfileName,"w")
intervalLength = int(args[2])
rateMatrixFilename=str(args[3])
groupMatrixFilename=str(args[4])
scale = 0
if len(args)==6:
	scale = int(args[5])
if len(args)==7:
	scalemin = int(args[5])
	scalemax=  int(args[6])
	scale = np.random.uniform(scalemin,scalemax)
	print "scale in (",scalemin,",",scalemax,"): ",scale
#np.random.seed(5)
outGroupMap = {}
inGroupMap = {}
def generateInterval(start, end, rateFileName, groupFileName):
	global groupMap
	rateFile = open(rateFileName,'r')
	groupFile = open(groupFileName,'r')
	ratematrix = []
	groupmatrix = []
	intLength = end-start
	vertexCount = 0
	for line in rateFile:
		rates = line.split()
		raterow = []
		for r in rates:	
			val = float(r)
			if scale>0:
				scalefact = (10.**float(scale/3.))/1000.
				val =scalefact*val	
				print "scalefact: ",scalefact," to ",val
			raterow.append(float(val))
		ratematrix.append(raterow)

	for line in groupFile:
		groups = line.split()
		grouprow = []
		for g in groups:
			grouprow.append(int(g))
		groupmatrix.append(grouprow)
	
	vertexIter = 0
	for i in range(0,len(groupmatrix)):
		for j in range(0,len(groupmatrix[0])):
			groupsize = groupmatrix[i][j]		
			vertexCount += groupsize
			for k in range(0,groupsize):
				outGroupMap[vertexIter] = i
				inGroupMap[vertexIter] = j
				vertexIter+=1
	
	for i in range(0,vertexCount):
		for j in range(0,vertexCount):
			#sample edges
			if i!=j:
				if ratematrix[outGroupMap[i]][inGroupMap[j]]>0:
					currentRate = ratematrix[outGroupMap[i]][inGroupMap[j]]
					nextTime = 2+np.random.exponential(1./currentRate)
					while nextTime<intLength:

						outTime=int(start+nextTime)
						outfile.write(str(i)+" "+str(j)+" "+str(outTime)+"\n")
						timeInt = np.random.exponential(1./currentRate)
						#print "time drawn: ",timeInt," from rate: ",1./currentRate," total: ",nextTime," adjusted: ",start+nextTime
						nextTime += timeInt
				

	rateFile.close()
	groupFile.close()



generateInterval(0,intervalLength,rateMatrixFilename,groupMatrixFilename)


