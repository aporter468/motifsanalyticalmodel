
#Usage: generateGraphsetsManualMultiWindow.sh ratefile groupfile
#same as generateGraphsetsManual.sh but repeats over multiple windows with different rate factors randomized by the givenGroupsSynthesizer script and combines them

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


outdir="manualGraphsetMultInt"$num
mkdir $outdir
#top of seq is # of graphs you want
for i in `seq 1 1`;
do
	subdir=$outdir"/set"$i
	mkdir $subdir
	fulloutgraph=$subdir"/fullGraph.txt"
	#top of seq is # of intervals on the graph each with length $windowLength
	for j in `seq 1 32`;
	do
		outgraph=$subdir"/syntheticGraph_"$j".txt"
		python givenGroupsSynthesizer.py givengroupraw.txt $windowLength $manratefile $mangroupfile 5 8
		sort -n -k3 givengroupraw.txt > givengroupsorted.txt
		python data_renamer.py givengroupsorted.txt > $outgraph
		tshift=$(($(($j-1))*$windowLength))
		echo "tshift: "$tshift
		cat $outgraph | awk -v t=$tshift '{print $1" "$2" "$3+t}' >> $fulloutgraph
	done

done
