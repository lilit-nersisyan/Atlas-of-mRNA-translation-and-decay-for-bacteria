library(Rsubread)

filedir <- "/crex/proj/snic2019-30-56/sample_files_mj/bsub/count_samples/bam/KD_batch2/"
setwd(filedir)
gr <- "/crex/proj/snic2019-30-56/mengjun/genome.annotation/bacteria/bsub.featureCount.txt"
files <- list.files(pattern=".bam$")
tmpcount <- c()
for (i in files){
  x <- featureCounts(i, isGTFAnnotationFile=F, annot.ext = gr, useMetaFeatures=F, strandSpecific = 1, 
                     allowMultiOverlap = T, minOverlap=1, countMultiMappingReads = F, isPairedEnd=F, 
                     requireBothEndsMapped=F, countChimericFragments = F, largestOverlap=T)
  
  
  tmpcount <- cbind(tmpcount,x$counts)
}

bsub.featureCount.KD.hs <- tmpcount
save(bsub.featureCount.KD.hs, file = "/crex/proj/snic2019-30-56/mengjun/microbiome_revision/reads_count/bsub.featureCount.KD.hs.rda")
