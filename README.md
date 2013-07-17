How to get started with BACON:
==============================

A scientific paper describing the algorithm and some results can be found in [PLoS ONE at this link](http://dx.plos.org/10.1371/journal.pone.0068358).


GNU Octave set-up
-----------------

Begin with a working installation of GNU Octave (which is generally Matlab-compatible) with following packages installed:

* general -- for the function parcellfun() for multi-core processing
* gsl -- for the function psi()
* miscellaneous
* optim
* statistics
* struct

The data
--------

The repo provides data from the [DREAM4 in silico challenge](http://wiki.c2b2.columbia.edu/dream/index.php/D4c2), which is the data set used in the PLoS ONE paper. Given these data in the ./data directory, the Octave script dream4script.m will infer interactions from the data in a number of configuations. As written, the script will use 2 processor cores to infer interactions on each network, first on all 5 time-series for the network and then on each time-series individually, with 10 random starts each. A final results object (.mat) is saved in the working directory. This file is large and has not been provided here.


Running the algorithm
---------------------

To obtain the "gold standard" networks for performance evaluation, run processGoldStandards.R in R (or, the results have already been calculated and are provided in the directory ./goldstdmats). These are the "gold standard" networks used to evaluate algorithm performance after inference, and are given in a different format in the DREAM4 data. See script for details.

The script dream4resultsCollection.m reads the results object and gold standard matrices from above and computes the AUROC and AUPR for each case, and further summarizes performance results. These results and summaries are normally written to the working directory as .mat files, and here they have been provided in the ./results10gene directory.

Finally, the code file postanalysis.m gives some basic comparisons between the clustering and no-clustering models, based on the .mat files written above. This is not a script, but provides all necessary final results in variable form.


Questions? 
----------

Please contact me through:
https://github.com/briangodsey/bacon-for-genetic-networks


Want to contribute?
-------------------

Feel free to fork this repo, submit pull requests, etc! I'm always happy to collaborate.



