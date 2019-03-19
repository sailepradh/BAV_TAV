if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")

library("edgeR")


#counttable= read.table( "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/SP4_p0.05/corrected_CNVs/Count_nlps_wlps_2.txt",
                        #sep ="\t", header= FALSE, row.names=1)


#setwd("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_callsFeb/Differential_Interaction_PP/")


setwd("//Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected/DL_PP/")
counttable_raw = read.table( "BAVTAV.Proximities_SP4_p01_filtered_corrected.txt",
                        sep ="\t",
                        header= FALSE
                        , row.names=1)

counttable <- counttable_raw [,c(4:9)]
head (counttable)
colnames(counttable) <- c("TAV_rep1", "TAV_rep2","TAV_rep3","BAV_rep1", "BAV_rep2", "BAV_rep3")

mobDataGroups <- c("TAV", "TAV", "TAV",
                   "BAV", "BAV", "BAV")
d <- DGEList(counts=counttable,group=factor(mobDataGroups))

dim(d)
d.full <- d
head(d$counts)

#head(cpm(d))
#apply(d$counts, 5, sum)


##### Very very  lenient filter applied in here dont expect to get a lot of loops after this filer
#keep <- rowSums(cpm(d)>5) >= 2
#d <- d[keep,]
#dim(d)21hicap_expr_PP_raw_18_Sep_2018.py

## Libarary size of the individual experiments
#TAV2431 = 2129053 + 104216145
#TAV2515 = 1975040 + 94842378
#TAV2709 = 9256828 + 83131288
#BAV2375 = 2762832 + 132399792
#BAV2424 = 2153215 + 115494826
#BAV2714 = 8315224 + 76183712


TAV2431 = 2909990 + 170232550
TAV2515 = 2934978 + 159740064
TAV2709 = 10677549 + 160319318
BAV2375 = 3945034 + 222762884
BAV2424 = 3151362 + 193090408
BAV2714 = 9605814 + 147196693
 


d$samples$lib.size <- c(TAV2431,TAV2515,TAV2709,BAV2375, BAV2424,BAV2714)
d$samples


d <- calcNormFactors(d)
d

par(mar=c(6, 6, 6, 6) + 0.8)
plotMDS(d, method ="bcv",
        col=as.numeric(d$samples$group),
        cex.lab =  1.5,
        ylim = c(-1.0, 1.0),
        main = "\n MDS plot of Significant Interaction in BAV and TAV Replicates\n",
        cex.main = 2.0
        )

legend("topright", sort(as.character(unique(d$samples$group))), col=1:2, pch=20, cex =1.5 )

d1 <- estimateCommonDisp(d, verbose=T)
d1$common.dispersion
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method ="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
d2
plotBCV(d2)

#################################################
## testing with exact model
et12 <- exactTest(d1, pair=c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.1)
summary(de1)

de1tags12 <- rownames(d1)[as.logical(de1)]
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")

#################################################
## Real stuff  on glm model 
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(-1,1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.1, lfc = 2)
#de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.1)
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
#xticks <- seq(-8,  8, 2)
yticks <- seq(-8, 8, 2)
axis(1)
axis(2 ,at = yticks, labels = yticks,  las=2)

abline(h = c(-2, 2), col = "blue")
legend("topright",  legend = c("TAV", "Both","BAV"),   fill= c("darkgreen","grey","darkmagenta"),
       cex = 1.5,
       bty = "n")

summary(de2)

DE_Loops <- lrt12$table
DE_Loops$p_adjust <- p.adjust(DE_Loops$PValue, method = "BH")
dim(DE_Loops)
dim(counttable_raw)
result_df <- merge.data.frame(counttable_raw , DE_Loops,
                              by="row.names", 
                              all.x=TRUE )

DE_Loop_significant <- subset(result_df, result_df$p_adjust < 0.1)

result_df <- merge.data.frame(counttable_raw , DE_Loops,
              by="row.names", 
              all.x=TRUE )


write.table(DE_Loop_significant,
            "DL_significant_03_06.tsv",
            sep="\t", quote = F, row.names = FALSE)

write.table(result_df,
            "DL_Loop_03_06.tsv",
            sep="\t", quote = F, row.names = FALSE)


