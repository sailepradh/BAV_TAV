library(edgeR)
library(statmod)
library(limma)
library(Rsubread)

setwd ("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/RNA_data")

#Read Counts and sample info
counts_1 <- read.csv("gene_count_matrix.csv",header=TRUE)
gene_2_ID <- read.csv("sorted_Gene_ID.csv",header=FALSE)
colnames (counts_1) <- c("Ensembl_ID",
                         "BAV2375_RNA_EC","BAV2424_RNA_EC","BAV2714_RNA_EC",
                         "TAV2431_RNA_EC","TAV2515_RNA_EC","TAV2709_RNA_EC")
colnames(gene_2_ID) <- c("Ensembl_ID", "Ensembl_ID_Genes_ID")
Both_combined <- merge (counts_1, gene_2_ID, by = "Ensembl_ID")
#write.table(Both_combined [,c(13,2:12)], "Counts_genes_ID.csv",sep=',',row.names = FALSE, quote = FALSE)


Counts <- Both_combined [,c(2,3,4,5,6,7)]
head(Counts)
row.names(Counts) <- Both_combined[,8]

DataGroups <- c ("BAV","BAV","BAV",
                 "TAV","TAV","TAV")

d <- DGEList(counts=Counts,group=factor(DataGroups))
d
dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>10) >= 3
d <- d[keep, keep.lib.sizes=FALSE]
dim(d)

AveLogCPM <- aveLogCPM(dgeObj, gene.length = tx_lens)
hist(AveLogCPM)

d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d
apply(d$counts, 2, sum)

plotMDS(d, top =500, col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12, n=10)

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)


de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")


design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1))
topTags(lrt12, n=10)

de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.1)
de2tags12 <- rownames(d2)[as.logical(de2)]
summary(de2)
de2tags12 <- rownames(d2) [as.logical(de2)]
de2tags12_BAV <- rownames(d2) [as.logical(de2==-1)]
de2tags12_TAV <- rownames(d2) [as.logical(de2== 1)]


plotSmear(lrt12, de.tags= de2tags12,
          pch = 19 , col = "grey", cex=0.5, smearWidth=0.5)

par(mar=c(8, 8, 8, 8) + 1)
plotMD(lrt12, status=de2, values=c(1,0,-1),pch = 19,cex=0.5,
       col=c("darkgreen","grey","darkmagenta"),
       main ="Differential Interaction in two conditions",
       cex.main = 2.0 ,
       cex.lab = 1.5, axes = FALSE,
       legend = FALSE)


yticks <- seq(-8, 8, 2)
axis(1)
axis(2 ,at = yticks, labels = yticks,  las=2)

abline(h = c(-2, 2), col = "blue")
legend("topright",  legend = c("TAV", "Both","BAV"),   fill= c("darkgreen","grey","darkmagenta"),
       cex = 1.5,
       bty = "n")

summary(de2)

DE <- lrt12$table
DE$p_adjust <- p.adjust(DE$PValue, method = "BH")
dim(DE)

write.table(DE,
            "Differential_expression.txt",
            sep="\t", quote = F, row.names = TRUE)
