# Exploration and retrieval of whole-metagenome sequencing samples

## Supporting files for generating simulated data using metaSim 
* use taxon profile given in bioinformaticsRevision.mprf
* use abundance values given in <experiment>-label-abundance files
* additional details available in paper

## Matlab code for reproducing figures
* run runPlots.m to reproduce figures in the paper
 * check for dependencies in configPaths.m 
* each <data>_<condition>.mat file contains 3 matrices for specific k-mer (K), unspecified k-mer (S), and FIGfam (F)
 * dimensions of recordK - (3 matrices), (k-mer 12 21 30), (entropies), (not needed), (query)
 * dimensions of recordS - (3 matrices), (entropies), (not needed), (query)
 * dimensions of recordF - (3 matrices), (entropies), (length 21 30), (query)
 * each figure*randomTest.mat file contains results of the random test
 * abundanceLabel contains abundance values estimated using MetaPhlAn
* More descriptions available within the code

## Comments
* synthetic data names: bioRev - HIGH-C, bioRev-2 - HIGH-VAR, bioRev-3 - LOW-C, bioRev-4 - MIXED-C
