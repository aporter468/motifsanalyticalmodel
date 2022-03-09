See runemail.sh for how to run a trial on the email data set
See Synthesizer/syntheticTest.sh for runtime comparison test (update script with path to Snap installation to run)

Datasets:
email_clean_500day.txt is EU research institute relabeled and edges in the first 500 days only selected.

email_clean_500day_thresholdclusternodes.txt is the list of nodes with degree >10% the maximum static degree in email_clean_500day.txt

email_thresholdcluster_relabeledgraph.txt is the (relabeled) graph induced by the largest cluster in the graph induced by email_clean_500day_thresholdclusternodes.txt.
