ratefile="timingtest_TASBMparams/rates.txt"
groupfile="timingtest_TASBMparams/groups.txt"
pathtosnap="../../Snap-4.1/examples/temporalmotifs"
groups=6
delta=1000
T=1000
graph="synthgraph.txt"
testtype=0
./main_gen $ratefile $groupfile $groups $delta $T $graph $testtype > testout.txt

if [ "$testtype" -eq "1" ];
then
#	sort -k 3,3 -n $graph >$graph"_sorted"
	$pathtosnap"/temporalmotifsmain" -i:$graph -o:"result.txt" -delta:$delta
fi
