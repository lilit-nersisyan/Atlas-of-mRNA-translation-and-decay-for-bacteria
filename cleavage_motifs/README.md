# Source data and script to generate cleavage motif logo files in Figure 1 

The 5´P cleavage motifs were computed from the sequence composition of transcripts, involving the region 4 nucleotides upstream and downstream from 5´ mapping positions of all the reads. If multiple reads map to the same position, we consider the base composition of that region multiple times as well. Using the obtained base frequencies, we then proceed to producing the sequence logos with the R package ggseqlogo, using Shannon entropy (bits) to compute the contribution of each nucleotide at each position.

The source data are the fivpeseq output files named "codon_pauses.txt" that contain 5' endpoint counts at a given distance from each codon. 
The script motif_codons.Rmd takes these files as input and generates logo files of cleavage motifs, taking the codonwise cumulative counts +/- 4 bases around the distance 0. The results are written to PDF files. 
