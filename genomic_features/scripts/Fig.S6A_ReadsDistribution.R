library(CAGEfightR)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tidyverse)
load("/data/bsub.hierachy.nc.rda")
load("/data/ecoli.hierachy.nc.rda")

cbPalette <- c("#999999", "#F0E442", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#000000", "#9999CC","#0072B2" )

GenomicDistribution_reads <- function(bw.dir, plus.pattern, minus.pattern,
                                      tx_model, outputprefix, tx_length, namebase){
  old.dir <- getwd()
  bw_plus.files <- list.files(bw.dir, pattern = "plus.bw")
  print(bw_plus.files)
  bw_minus.files <- list.files(bw.dir, pattern = "minus.bw")
  res <- c()
  setwd(bw.dir)
  for(i in namebase){
    print(i)
    print(bw_plus.files[grep(i, bw_plus.files)])
    plus.tmp <- import(bw_plus.files[grep(i, bw_plus.files)])
    strand(plus.tmp) <- "+"
    minus.tmp <- import(bw_minus.files[grep(i, bw_minus.files)])
    strand(minus.tmp) <- "-"
    bw.tmp <- c(plus.tmp, minus.tmp)
    trial <- assignTxType(bw.tmp, tx_model)
    trial2 <- data.frame(trial)
    trial2 <- group_by(trial2, by=txType)
    trial2 <- data.frame(mutate(trial2, totalSig=sum(score)))
    trial2.sub <- trial2[,c("txType","totalSig")]
    trial2.sub <- trial2.sub[!duplicated(trial2.sub),]
    trial2.sub <- trial2.sub[trial2.sub$txType!="rRNA", ]
    res.tmp <- data.frame(region=trial2.sub$txType,
                          count=trial2.sub$totalSig, 
                          percenctage=as.vector(trial2.sub$totalSig)/sum(as.vector(trial2.sub$totalSig))*100)
    
    res.tmp$group <- i
    res <- rbind(res, res.tmp)
  }
  genomicLevels <- c("TSS", "UTR5","cds","UTR3","UTRinternal","tRNA", "otherRNA", "intergenic")
  genomic.tmp <- unique(res$region)
  genomicLevels <- genomicLevels[genomicLevels%in% genomic.tmp]	
  res$region <- factor(res$region, levels = genomicLevels)
  res$group <- factor(res$group, levels = namebase)
  p <- ggplot(data=res, aes(x=group, y=percenctage, fill=region)) +
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = cbPalette)
  
  pdf(file=paste0(outputprefix, "Reads_distribution_byRegions.pdf"))
  print(p)
  dev.off()
  setwd(old.dir)
  res
  
}


bw.dir <- "/crex/proj/snic2019-30-56/sample_files_mj/bsub/pooled_samples/"
plus.pattern="pooled_plus.bw"
minus.pattern="pooled_minus.bw"
namebase <- c("ctrl","rnjA_","rnjB_","rny_", "^hs_","mup",
              "cam_","cam-ms","cam-ps","rnjA-hs","rnjB-hs","rny-hs",
              "salt","stat-24","stat-48","stat-d8")

##
outputprefix <- "/crex/proj/snic2019-30-56/mengjun/microbiome_revision/reads_distribution/bsub_"

reads.dist.bsub <- GenomicDistribution_reads(bw.dir=bw.dir, plus.pattern=plus.pattern, 
                                             minus.pattern=minus.pattern,
                                             tx_model=bsub.hierachy.nc, 
                                             outputprefix=outputprefix, 
                                             tx_length=bsub.genomicRegion.width, 
                                             namebase=namebase)


bw.dir <- "/crex/proj/snic2019-30-56/sample_files_mj/ecol/pooled_samples/"
plus.pattern="pooled_plus.bw"
minus.pattern="pooled_minus.bw"
namebase=c("^ctr-dh5a","rnja-ctr-dh5a","^cam-dh5a", "rnja-cam-dh5a",
           "ctr-mg1655","hs-mg1655","mup-mg1655","cam-mg1655")

##
outputprefix <- "/crex/proj/snic2019-30-56/mengjun/microbiome_revision/reads_distribution/ecoli_hierachy_"
reads.dist.ecoli <- GenomicDistribution_reads(bw.dir=bw.dir, plus.pattern=plus.pattern, 
                                              minus.pattern=minus.pattern,
                                              tx_model=ecoli.hierachy.nc, 
                                              outputprefix=outputprefix, 
                                              tx_length=bsub.genomicRegion.width, 
                                              namebase)
