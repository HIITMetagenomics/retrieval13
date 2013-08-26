# retrieval13



## Matlab code for reproducing figures
To generate samples from synthetic data, use generateSamples.m
* Install metaSim and Jellyfish, check for dependencies explained in the code
* .mprf files for LOW and HIGH are avaliable
* Abundance profile of each sample and their class label is available in [LOW/HIGH].abundance.label
    
To reproduce the figures in the paper run fig*.m
* Check for dependencies in configPaths.m 
* Each <data>.mat file contains 3 matrices for specific k-mer (K), unspecified k-mer (S), and FIGfam (F)
 * Dimensions of recordK - (3 metrices), (k-mer 12 21 30), (entropies), (not needed), (query)
 * Dimensions of recordS - (3 metrices), (entropies), (not needed), (query)
 * Dimensions of recordF - (3 metrices), (entropies), (length 21 30), (query)
 * Each figure*randomTest.mat file contains results of the random test
  * Description is available within the code
