Synthesizers:

- All methods use givenGroupSynthesizer.py
	usage (windowLength is integer, rate, group files square matrices of same dimensions):
	 python givenGroupsSynthesizer.py <tempfile name>  <window length> <rate file>  <group file>
	with optional arguments:
	<scale> or <scale min> <scale max> 

	
- Batch generating scripts use givenGroupsSynthesizer to generate sets of random graphs and incorporate graph-level variables, i.e. changing over time or scaling rates
	usage (set other variables in scripts):
	<batchscript choice> <ratefile> <groupfile> 
