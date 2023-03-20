library(edgeR)

#######################rna-kd.hs#####################
bsub.anno <- read.table("/data/bsub.anno.txt", header=T)
load("/data/bsub.featureCount.KD.hs.rda")
replicates <- c(paste0("rep",1:3), paste0("rep",1:3), 
                paste0("rep",1:3), paste0("rep",1:3),
                paste0("rep",1:3),paste0("rep",1:3),
                paste0("rep",1:3))
groups <- c(rep("rnja",3), rep("rnja-hs",3), 
            rep("rnjb",3), rep("rnjb-hs",3),
            rep("rny",3),rep("rny-hs",3),
            rep("ctrl",3))
sample.names <- colnames(bsub.featureCount.KD.hs)
library.size <- colSums(bsub.featureCount.KD.hs)

RNAseq_design <- data.frame(replicates, groups, sample.names, 
                            library.size = library.size)
y <- DGEList(counts=bsub.featureCount.KD.hs, samples = RNAseq_design, 
             lib.size = RNAseq_design$library.size)

keep <- rowSums(y$counts> 1) >= 3
y <- y[keep, , keep.lib.sizes=T]

groups2 <- factor(y$sample$group)
design <- model.matrix(~0+groups2, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)


##quasi-likelihood: better control of type I error
fit2 <- glmQLFit(y, design)
DE_result_rnja <- glmQLFTest(fit2, contrast = c(-1,1,0,0,0,0,0))
DE_result_rnja <- topTags(DE_result_rnja, n=nrow(DE_result_rnja$table))
DE_result_rnja <- DE_result_rnja$table[,c("logFC","FDR")]
colnames(DE_result_rnja) <- paste0("rnja_",colnames(DE_result_rnja))

DE_result_rnja.hs <- glmQLFTest(fit2, contrast = c(-1,0,1,0,0,0,0))
DE_result_rnja.hs <- topTags(DE_result_rnja.hs, n=nrow(DE_result_rnja.hs$table))
DE_result_rnja.hs <- DE_result_rnja.hs$table[row.names(DE_result_rnja),c("logFC","FDR")]
colnames(DE_result_rnja.hs) <- paste0("rnja.hs_",colnames(DE_result_rnja.hs))

DE_result_rnjb <- glmQLFTest(fit2, contrast = c(-1,0,0,1,0,0,0))
DE_result_rnjb <- topTags(DE_result_rnjb, n=nrow(DE_result_rnjb$table))
DE_result_rnjb <- DE_result_rnjb$table[row.names(DE_result_rnja),c("logFC","FDR")]
colnames(DE_result_rnjb) <- paste0("rnjb_",colnames(DE_result_rnjb))

DE_result_rnjb.hs <- glmQLFTest(fit2, contrast = c(-1,0,0,0,1,0,0))
DE_result_rnjb.hs <- topTags(DE_result_rnjb.hs, n=nrow(DE_result_rnjb.hs$table))
DE_result_rnjb.hs <- DE_result_rnjb.hs$table[row.names(DE_result_rnja),c("logFC","FDR")]
colnames(DE_result_rnjb.hs) <- paste0("rnjb.hs_",colnames(DE_result_rnjb.hs))

DE_result_rny <- glmQLFTest(fit2, contrast = c(-1,0,0,0,0,1,0))
DE_result_rny <- topTags(DE_result_rny, n=nrow(DE_result_rny$table))
DE_result_rny <- DE_result_rny$table[row.names(DE_result_rnja),c("logFC","FDR")]
colnames(DE_result_rny) <- paste0("rny_",colnames(DE_result_rny))

DE_result_rny.hs <- glmQLFTest(fit2, contrast = c(-1,0,0,0,0,0,1))
DE_result_rny.hs <- topTags(DE_result_rny.hs, n=nrow(DE_result_rny.hs$table))
DE_result_rny.hs <- DE_result_rny.hs$table[row.names(DE_result_rnja),c("logFC","FDR")]
colnames(DE_result_rny.hs) <- paste0("rny.hs_",colnames(DE_result_rny.hs))

res.diff <- cbind(DE_result_rnja, DE_result_rnjb, DE_result_rny,
                  DE_result_rnja.hs, DE_result_rnjb.hs, DE_result_rny.hs)

res.diff.kd.only <- cbind(DE_result_rnja, DE_result_rnjb, DE_result_rny)

######sig
DE_result_rnja.sig <- DE_result_rnja[DE_result_rnja$rnja_FDR < 0.05, ]
DE_result_rnjb.sig <- DE_result_rnjb[DE_result_rnjb$rnjb_FDR < 0.05, ]
DE_result_rny.sig <- DE_result_rny[DE_result_rny$rny_FDR < 0.05, ]
DE_result_rnja.hs.sig <- DE_result_rnja.hs[DE_result_rnja.hs$rnja.hs_FDR < 0.05, ]
DE_result_rnjb.hs.sig <- DE_result_rnjb.hs[DE_result_rnjb.hs$rnjb.hs_FDR < 0.05, ]
DE_result_rny.hs.sig <- DE_result_rny.hs[DE_result_rny.hs$rny.hs_FDR < 0.05, ]

all.sig <- unique(c(row.names(DE_result_rnja.sig),row.names(DE_result_rnjb.sig),
                    row.names(DE_result_rny.sig),row.names(DE_result_rnja.hs.sig),
                    row.names(DE_result_rnjb.hs.sig),
                    row.names(DE_result_rny.hs.sig)))
kd.sig <- unique(c(row.names(DE_result_rnja.sig),row.names(DE_result_rnjb.sig),
                   row.names(DE_result_rny.sig)))

res.diff.sig <- res.diff[row.names(res.diff) %in% all.sig,]
res.diff.KD.only.sig <- res.diff.kd.only[row.names(res.diff.kd.only) %in% kd.sig,]
res.diff.sig.FC <- res.diff.sig[,c(1,3,5,7,9,11)]
res.diff.KD.only.sig.FC <- res.diff.KD.only.sig[,c(1,3,5)]

res.diff.sig.FDR <- res.diff.sig[,c(2,4,5,8,10,12)]
res.diff.sig$id <- row.names(res.diff.sig)
bsub.anno <- bsub.anno[bsub.anno$genomic.feature!="rRNA",]
res.diff.sig <- merge(res.diff.sig, bsub.anno)
res.diff.sig <- res.diff.sig[order(res.diff.sig$id),]
write.csv2(res.diff.sig, file="rnaseKD.hs.significant.changed.csv")

res.diff.KD.only.sig$id <- row.names(res.diff.KD.only.sig)
res.diff.KD.only.sig <- merge(res.diff.KD.only.sig, bsub.anno)
res.diff.KD.only.sig <- res.diff.KD.only.sig[order(res.diff.KD.only.sig$id),]
write.csv2(res.diff.KD.only.sig, file="rnaseKD.only.significant.changed.csv")

###
No.sig.genomicFeatures <- data.frame(table(res.diff.sig$genomic.feature))
No.all.genomicFeatures <- data.frame(table(bsub.anno$genomic.feature))
No.regulated.genomicFeatures <- No.all.genomicFeatures
colnames(No.regulated.genomicFeatures) <- c("genomic.feature", "all")
No.regulated.genomicFeatures$significant.changed <- No.sig.genomicFeatures$Freq
write.table(No.regulated.genomicFeatures, file="No.regulated.genomicFeatures.KD.hs.txt",quote=F,
            sep="\t", row.names=F)

###
No.sig.KD.only.genomicFeatures <- data.frame(table(res.diff.KD.only.sig$genomic.feature))
No.all.genomicFeatures <- data.frame(table(bsub.anno$genomic.feature))
No.regulated.genomicFeatures <- No.all.genomicFeatures
colnames(No.regulated.genomicFeatures) <- c("genomic.feature", "all")
No.regulated.genomicFeatures$significant.changed <- No.sig.KD.only.genomicFeatures$Freq
write.table(No.regulated.genomicFeatures, file="No.regulated.genomicFeatures.KD.only.txt",quote=F,
            sep="\t", row.names=F)

library("RColorBrewer")
res.diff.sig.FC2 <- round(res.diff.sig.FC[row.names(res.diff.sig.FC) %in%res.diff.sig$id,])
res.diff.KD.only.sig.FC2 <- round(res.diff.KD.only.sig.FC[row.names(res.diff.KD.only.sig.FC) %in% res.diff.KD.only.sig$id,])


paletteLength <- 9
#colors (one less than breaks
colfunc <- colorRampPalette(c("navyblue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(res.diff.sig.FC2), 0, length.out=ceiling(paletteLength/2)), 
              seq(max(res.diff.sig.FC2)/paletteLength, max(res.diff.sig.FC2), length.out=floor(paletteLength/2)))

pheatmap(res.diff.sig.FC2, color=colfunc,cluster_rows = T,
         cluster_cols = T, show_rownames = F, breaks=myBreaks,
         clustering_method="ward.D2")

utr5 <- bsub.anno[bsub.anno$genomic.feature=="5utr",]$id
smallRNA <- bsub.anno[bsub.anno$genomic.feature=="smallRNA",]$id

###5utr
res.diff.sig.FC2.utr5 <- res.diff.sig.FC2[row.names(res.diff.sig.FC2)%in%utr5,]
pheatmap(res.diff.sig.FC2.utr5, color=colfunc,cluster_rows = T,
         cluster_cols = T, show_rownames = F, breaks=myBreaks,
         clustering_method="ward.D2", filename="bsub_mutants_hs_FC_heatmap_utr5.pdf")

myBreaks <- c(seq(min(res.diff.KD.only.sig.FC2), 0, length.out=ceiling(paletteLength/2)), 
              seq(max(res.diff.KD.only.sig.FC2)/paletteLength, max(res.diff.KD.only.sig.FC2), length.out=floor(paletteLength/2)))
res.diff.KD.only.sig.FC2.utr5 <- res.diff.KD.only.sig.FC2[row.names(res.diff.KD.only.sig.FC2)%in%utr5,]

pheatmap(res.diff.KD.only.sig.FC2.utr5, color=colfunc,cluster_rows = T,
         cluster_cols = T, show_rownames = F, breaks=myBreaks,
         clustering_method="ward.D2", filename="bsub_mutants_only_FC_heatmap_utr5.pdf")

###smallRNA
res.diff.sig.FC2.smallRNA <- res.diff.sig.FC2[row.names(res.diff.sig.FC2)%in%smallRNA,]
pheatmap(res.diff.sig.FC2.smallRNA, color=colfunc,cluster_rows = T,
         cluster_cols = T, show_rownames = F, breaks=myBreaks,
         clustering_method="ward.D2", filename="bsub_mutants_hs_FC_heatmap_miscRNA.pdf")

res.diff.KD.only.sig.FC2.utr5 <- res.diff.KD.only.sig.FC2[row.names(res.diff.KD.only.sig.FC2)%in%smallRNA,]

pheatmap(res.diff.KD.only.sig.FC2.utr5, color=colfunc,cluster_rows = T,
         cluster_cols = T, show_rownames = F, breaks=myBreaks,
         clustering_method="ward.D2", filename="bsub_mutants_only_FC_heatmap_miscRNA.pdf")


####pca

data.pca <- cpm(y,normalized.lib.sizes = F, log=T)
cds <- bsub.anno[bsub.anno$genomic.feature=="cds",]$id
data.pca <- data.pca[row.names(data.pca) %in% cds,]
CountTable.filter.cpm.PC<-prcomp(t(data.pca), 
                                 center = TRUE, scale. = T)
variance_explained_all <- summary(CountTable.filter.cpm.PC)
libname1 <- unlist(lapply(row.names(CountTable.filter.cpm.PC$x), function(x)strsplit(x, split="_")[[1]][1]))
libname2 <- unlist(lapply(row.names(CountTable.filter.cpm.PC$x), function(x)strsplit(x, split="_")[[1]][2]))
libname <- paste0(libname1, ".", libname2)
PCi_all<-data.frame(CountTable.filter.cpm.PC$x, 
                    Library=libname)
ve <- round(variance_explained_all$importance[2,] *100,digits = 1)
p <- ggplot(PCi_all,aes(x=PC1,y=PC2,col=Library)) + geom_point(size=5) + xlab(paste0("PC1 (", ve[1],"%)")) + ylab(paste0("PC2 (", ve[2],"%)"))

pdf(file = "PCA_plot_mutants.hs.pdf", width=6, height=4)
print(p)
dev.off()

data.pca2 <- cpm(y,normalized.lib.sizes = F, log=F)
data.pca2 <- data.pca2[row.names(data.pca2) %in% cds,]
mutant.only <- data.pca2[,c(1:3,7:9,13:15,19:21)]
keep <- rowSums(mutant.only> 1) >= 3
mutant.only <- mutant.only[keep, ]

CountTable.filter.cpm.PC<-prcomp(t(log(mutant.only+1)), 
                                 center = TRUE, scale. = T)
variance_explained_all <- summary(CountTable.filter.cpm.PC)
libname1 <- unlist(lapply(row.names(CountTable.filter.cpm.PC$x), function(x)strsplit(x, split="_")[[1]][1]))
libname2 <- unlist(lapply(row.names(CountTable.filter.cpm.PC$x), function(x)strsplit(x, split="_")[[1]][2]))
libname <- paste0(libname1, ".", libname2)
PCi_all<-data.frame(CountTable.filter.cpm.PC$x, 
                    Library=libname)
ve <- round(variance_explained_all$importance[2,] *100,digits = 1)
p <- ggplot(PCi_all,aes(x=PC1,y=PC2,col=Library)) + geom_point(size=5) + xlab(paste0("PC1 (", ve[1],"%)")) + ylab(paste0("PC2 (", ve[2],"%)"))
print(p)

pdf(file = "PCA_plot_mutants.only.pdf", width=6, height=4)
print(p)
dev.off()
