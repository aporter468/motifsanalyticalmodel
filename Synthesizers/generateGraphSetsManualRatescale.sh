#Usage: generateGraphsetsManualRatescale.sh ratefile groupfile
#same as generateGraphsetsManual.sh but outer loop determines randomly chosen scale on the ratefile (param $i passed to givenGroupSynthesizer used in randomization)
windowLength=10000
windowCount=1
vertexCount=100
delta=1000
outfileInit="tempEdges.txt"
outfileSorted="tempEdges2.txt"
outfile="syntheticGraph.txt"
motifOutfile="syntheticGraph_motifs.txt"
predOutputFile="predictionoutput.txt"

#num=1
#manratefile="manualRates"$num".txt"
#mangroupfile="manualGroups"$num".txt"
manratefile=$1
mangroupfile=$2

outdir="manualGraphsetRate"$num
mkdir $outdir
for i in `seq 1 10`;
do
	subdir=$outdir"/set"$i
	mkdir $subdir
	for j in `seq 1 30`;
	do
		outgraph=$subdir"/syntheticGraph_"$j".txt"
		python givenGroupsSynthesizer.py givengroupraw.txt $windowLength $manratefile $mangroupfile $i
		sort -n -k3 givengroupraw.txt > givengroupsorted.txt
		python data_renamer.py givengroupsorted.txt > $outgraph
	done

done
