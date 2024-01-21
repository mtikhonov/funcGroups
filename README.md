# funcGroups
MATLAB code to reproduce figures in Zhao, Cordero &amp; Tikhonov (2024), "Linear-regression-based algorithms can succeed at identifying microbial functional groups despite the nonlinearity of ecological function"

Code by Yuanchen Zhao

To use, run RunMe with no argument. The code will generate all data-dependent figures in the paper.

The precomputed data files (includeed) will be used, if available. Remove (or rename) a data file to recompute it from scratch.

Code uses "parfor" and will use the parallel computing toolbox with the default local configuration, if available.
