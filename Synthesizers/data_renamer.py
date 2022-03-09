#renames nodes in order of appearance for conveniece/efficiency in C++ model implementation

import numpy as np
import matplotlib.pyplot as plt
import sys

def readAndSortEdges(file):
    edgeDict = {}
    reverseDict = {}
    newName = {}
    namecounter=0
    with open (file, "r") as myfile:
        for line in myfile:
            lineparts = line.split()
            edge = (lineparts[0],lineparts[1])
            time = float(lineparts[2])
	    u = edge[0]
	    v = edge[1]
	    if u not in newName.keys():
		newName[u] = namecounter
		namecounter+=1
  	    if v not in newName.keys():
		newName[v] = namecounter
		namecounter+=1
	    print newName[u],newName[v],int(time)
            #sortedEdges.append(edgeDict[key])


args = sys.argv
readAndSortEdges(str(args[1]))

