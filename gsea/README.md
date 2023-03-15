# Scripts and source files to produce the GSEA results presented in Supplementary Table 3

Functional analysis was performed on the lists of genes with differential abundance or with differential frame protection between stress and control conditions. The FPI for each transcript was computed as described above. The fold change of the FPIs is computed as the difference of mean FPIs between the stress/mutant condition and the untreated/wild type control for each transcript. Differential abundance of 5’ endpoint counts for each transcript was computed as the expression log2 fold change with the DESeq2 R package using adaptive Student's t prior shrinkage estimator (apeglm). 

Gene set enrichment analysis (GSEA) was performed with the R package WebGestaltR based on the fold change differential expression or FPI values. GMT formatted files for each species were obtained by modifying Uniprot protein annotations and used as enrichment databases. The annotation GFF files for each species were taken as the basis to generate the reference gene lists. The significance p values are computed with a hypergeometic test and false discovery rate (FDR) is used for multiple test correction. 

The input data files to computed expression and FPI fold changes are in the frame_counts.zip and tr_descriptors.zip. The reference files for GSEA in the go folder. 
The stress_FPI-DEG.Rmd is the script for producing the GSEA output. 
