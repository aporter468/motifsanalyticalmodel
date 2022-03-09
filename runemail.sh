#example: graph that goes up to time 10k with 62023 edges; use delta 1000, T = 0.1*10000 and do window of second tenth of edges 
#timelimit=43200000
#edgecount=307869
edgecount=43035
timelimit=43186360
numnodegroups=8
delta=4320000
startwindow=2 # in range 0....1/windowfrac-1
endwindow=$((startwindow+1))
windowfrac=0.1
starttime=23175
timeinterval=$((timelimit-starttime))
#ingraph="datasets/email_clean_500day.txt"
ingraph="datasets/email_thresholdcluster_relabeledgraph.txt"
outedges="email_edges_used.txt"
outmotifs="email_motifs.txt"

./main $timeinterval $edgecount $ingraph $numnodegroups $outedges $outmotifs $delta $startwindow $endwindow $windowfrac $starttime
