library(edgeR)
library(statmod)
library(limma)
library(Rsubread)

#Read Tx Lenghts
tx.lens <- as.matrix(read.csv("tx_lens.txt", row.names = "transcript_id", sep = "\t"))

#Read Counts and sample info
counts<-as.matrix(read.csv("transcript_count_matrix_day_3_7_rep2_rep3.csv",row.names="transcript_id"))
sampleinfo <- read.csv("sample_key_day_3_7_rep2_rep3.txt", sep = '\t', row.names = 1)

#Group Samples
group <- paste(sampleinfo$condition)
group <- factor(group)
table(group)
timepoint <- sampleinfo$condition
design <- model.matrix(~0 + timepoint)
design

#Create the dgeObj
dgeObj <- DGEList(counts, group = group)

#Filter those with low counts and expression
myCPM <- cpm(dgeObj, gene.length = tx.lens )
thresh <- myCPM > 0.5
keep <- rowSums(thresh) > 3
table(keep)
dgeObj <- dgeObj[keep, keep.lib.sizes=FALSE]

#Average Log CPM Histogram
AveLogCPM <- aveLogCPM(dgeObj, gene.length = tx_lens)
hist(AveLogCPM)

#Calculate Normalisation Factors
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

#Visualise Samples
plotMDS(dgeObj, labels=group, cex=0.75)
plotMD(dgeObj,column = 1)
abline(h=0, col="red", lty=2, lwd=2)

#Estimate Dispersion
dgeObj <- estimateDisp(dgeObj,design, robust=TRUE)
plotBCV(dgeObj)
fit <- glmQLFit(dgeObj, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)

#Differential Expression
diffex <- makeContrasts(timepointday3-timepointday7, levels=design)
res <- glmQLFTest(fit, contrast = diffex)
res2 <- topTags(res, n = nrow(res$table))
is.de <- decideTestsDGE(res,adjust.method = "none")
summary(is.de)

jpeg('MD_plot_day3_vs_day7_rep2_rep3.jpg')
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
dev.off()
write.table(res2$table,"resSig_day_3_7_rep2_rep3.txt",sep='\t')


#tr <- glmTreat(fit, contrast=diffex, lfc=log2(1.5))
#tr_table <- topTags(tr, n = 100)
#is.de <- decideTestsDGE(tr)
#summary(is.de)

