source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library("edgeR")


#counttable= read.table( "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/SP4_p0.05/corrected_CNVs/Count_nlps_wlps_2.txt",
                        #sep ="\t", header= FALSE, row.names=1)



counttable= read.table( "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/DE_loop_test2/count_test.txt",
                        sep ="\t", header= FALSE, row.names=1)

colnames(counttable) <- c("nlps_rep1", "nlps_rep3","wlps_rep1", "wlps_rep3")
mobDataGroups <- c("LPS_Untreated", "LPS_Untreated", "LPS_Treated", "LPS_Treated")
d <- DGEList(counts=counttable,group=factor(mobDataGroups))

dim(d)
d.full <- d
head(d$counts)

#head(cpm(d))
#apply(d$counts, 2, sum)


##### Very very  lenient filter applied in here dont expect to get a lot of loops after this filer
#keep <- rowSums(cpm(d)>5) >= 2
#d <- d[keep,]
#dim(d)21hicap_expr_PP_raw_18_Sep_2018.py


## no of read pairs in corresponding
THP1nLPS1 = 2280581 + 79174043
THP1nLPS3 = 2155432 + 97124480
THP1wLPS1 = 3888185 + 126873916
THP1wLPS3 = 3613702 + 157841817

d$samples$lib.size <- c(81454624,99279912 , 130762101, 161455519)
d$samples


d <- calcNormFactors(d)
d

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("topright", as.character(unique(d$samples$group)), col=1:3, pch=20)

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

et12 <- exactTest(d1, pair=c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
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
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
#de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.1)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")
summary(de2)

DE_Loops <- lrt12$table
DE_Loops$p_adjust <- p.adjust(DE_Loops$PValue, method = "BH")
DE_Loop_significant <- subset(DE_Loops, DE_Loops$p_adjust < 0.05 )

write.table(DE_Loop_significant,
            "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/DE_loop_test2/DE_Loop_significant.tsv",
            sep="\t", quote = F, row.names = TRUE)

write.table(DE_Loops,
            "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/DE_loop_test2/DE_Loop_all.tsv",
            sep="\t", quote = F, row.names = TRUE)



#write.table(DE_Loop_significant,
#            "/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/SP4_p0.05/corrected_CNVs/DE_Loop_significant_2.tsv",
#            sep="\t", quote = F, row.names = TRUE)



#write.table(DE_Loops,
            #"/Volumes/Work_drive/prj/THP1_July_2018/data/raw_external/Sep_2018/SP4_p0.05/corrected_CNVs/DE_Loops.tsv",
            #sep="\t", quote = F, row.names = TRUE)

