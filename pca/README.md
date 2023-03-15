To generate a "ribosome protection phenotype" we took the sum of counts positioned 30 to 1 nucleotide upstream of each amino acid and concatenated the per amino acid scaled counts to obtain a vector that describes ribosome protection in each sample. These vectors were used as input for principal component analysis (PCA) performed with the prcomp function of the R package stats (v 3.6.1). The PCA plots were generated with the autoplotly package (v0.1.2) in R. 

The source data for the PCA calculations is an R data file *master.mat.codon-v2.RData* that has been generated based on fivepseq output files. It contains the following data columns: 

Taxonomy: "Genus"     "Phylum"     "taxon"      
Sample info: "sample"     "trt"        "trt_simple" "rep"        "source"   "ps"    "libsize"    
Fivepseq stats: "fft_score"  "frame_pref" "pattern"    "F0"         "F1"         "F2"   
5' counts at given distance from each codon: "ALA_GCA-30" "ALA_GCA-29" "ALA_GCA-28" "ALA_GCA-27" "ALA_GCA-26" "ALA_GCA-25" "ALA_GCA-24" "ALA_GCA-23" *etc.*

The script *source_script-pca.Rmd* takes this master data file as input, filters and subsets it for each figure and uses the *plot.pca* function based on ggfortify and autoplotly to produce the PCA plots. 
