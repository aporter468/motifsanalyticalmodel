#Usage: generateGraphsetsManual.sh ratefile groupfile
#Generates set of synthetic graphs from manual groups + manual rates file
#in rates + group files row is "out" group or "from" group, col is "in" or "to" group
windowLength=10000
windowCount=1
vertexCount=720
delta=1000
outfileInit="tempEdges.txt"
outfileSorted="tempEdges2.txt"
outfile="syntheticGraph.txt"
motifOutfile="syntheticGraph_motifs.txt"
predOutputFile="predictionoutput.txt"

#alt system to inline: set num and using below convention for rates/groups files
#num=1
#manratefile="manualRates"$num".txt"
#mangroupfile="manualGroups"$num".txt"

manratefile=$1
mangroupfile=$2
outdir="manualGraphset"
mkdir $outdir

#change top of seq to # of graphs you want
for j in `seq 1 5`;
do
	outgraph=$outdir"/syntheticGraph_"$j".txt"
	python givenGroupsSynthesizer.py givengroupraw.txt $windowLength $manratefile $mangroupfile
	sort -n -k3 givengroupraw.txt > givengroupsorted.txt
	python data_renamer.py givengroupsorted.txt > $outgraph
done


