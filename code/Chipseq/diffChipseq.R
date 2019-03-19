## DiffBind are differentially bound between two or more groups
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind", version = "3.8")

library (DiffBind)
tamoxifen <- dba(sampleSheet="tamoxifen.csv")
tamoxifen <- dba.count(tamoxifen)
tamoxifen <- dba.contrast(tamoxifen)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen.DB <- dba.report(tamoxifen)

samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "tamoxifen.csv"))
names(samples)
samples

tamoxifen <- dba(sampleSheet="tamoxifen.csv", dir=system.file("extra", package="DiffBind"))
t1 <- data(tamoxifen_peaks)

dba.plotHeatmap(tamoxifen)

data(tamoxifen_counts)
tamoxifen <- dba.count(tamoxifen, summits=250)
tamoxifen

plot(tamoxifen)

tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen

plot(tamoxifen, contrast=1)
tamoxifen.DB <- dba.report(tamoxifen)
tamoxifen.DB 


data(tamoxifen_analysis)
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)


dba.plotPCA(tamoxifen, contrast=1,label=DBA_TISSUE)

dba.plotMA(tamoxifen)
dba.plotVolcano(tamoxifen)


sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)

pvals <- dba.plotBox(tamoxifen)
corvals <- dba.plotHeatmap(tamoxifen)

data (tamoxifen.DB)
corvals <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE)


data (tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen,categories=DBA_CONDITION, block=tamoxifen$masks$MCF7)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen


dba.plotMA(tamoxifen,method=DBA_DESEQ2_BLOCK)


tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=T)[9:12]

tam.block <- dba.report(tamoxifen,method=DBA_ALL_METHODS_BLOCK,bDB=TRUE,bAll=TRUE)
tam.block
dba.plotVenn(tam.block,1:4,label1="edgeR",label2="DESeq2",
             label3="edgeR Blocked", label4="DESeq2 Blocked")

data(tamoxifen_peaks)
olap.rate <- dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

names(tamoxifen$masks)
dba.overlap(tamoxifen,tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive,
            mode=DBA_OLAP_RATE)
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)

tamoxifen_consensus <- dba.peakset(tamoxifen, consensus=c(DBA_TISSUE,DBA_CONDITION),
                                   minOverlap=0.66)


tamoxifen_consensus <- dba.peakset(tamoxifen, consensus=c(DBA_TISSUE,DBA_CONDITION),
                                    minOverlap=0.66)
tamoxifen_consensus <- dba(tamoxifen_consensus, mask=tamoxifen_consensus$masks$Consensus,
                           minOverlap=1)

tamoxifen_consensus
consensus_peaks <- dba.peakset(tamoxifen_consensus, bRetrieve=TRUE)
data(tamoxifen_counts)
tamoxifen <- dba.count(tamoxifen, peaks=consensus_peaks)
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_TISSUE, minOverlap=0.66)
dba.plotVenn(tamoxifen, tamoxifen$masks$Consensus)

data(tamoxifen_peaks)
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_CONDITION, minOverlap=0.33)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)

tamoxifen.OL <- dba.overlap(tamoxifen, tamoxifen$masks$Consensus)
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB

tamoxifen <- dba.peakset(tamoxifen,tamoxifen$masks$Consensus,
                        minOverlap=1,sampID="OL Consensus")
tamoxifen <- dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,
                          minOverlap=3,sampID="Consensus_3")
dba.plotVenn(tamoxifen,14:15)
data(tamoxifen_analysis)

tamoxifen.rep <- dba.report(tamoxifen,bCalled=TRUE,th=1)
onlyResistant <- tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant )
onlyResponsive <- tamoxifen.rep$Called2>=3 & tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups <- tamoxifen.rep$Called1>= 2 & tamoxifen.rep$Called2>=3
sum(bothGroups)
