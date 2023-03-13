# RNA-translation-and-decay-mapped-in-96-bacterial-species
Scripts used in the paper by Huch, Nersisyan et al. Nat. Micr. 2023 RNA  translation and decay mapped in 96 bacterial species

###	ribosome_protection_pca.Rmd and ribosome.master.mat.RData
Script to produce the PCA plots and the data containing one-dimensional vectors of amino-acid relative counts for each sample. The scripts loads the dataset when it is placed in the same directory and produces the plots.

###	gene_specific_frame_test.py
Script to perform statistical tests on the difference between frame-relative count distributions in each gene across samples. 

### graphlan_annotation.R and graphlan.data.RData 
Script produces text files with annotations that are used as input for the graphlan program to produce phylogenetic trees. The data contains values for periodicity, coverage, frame preference and identified enzymes and is loaded by the script if placed in the same directory. 

