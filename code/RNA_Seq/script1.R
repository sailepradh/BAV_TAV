library(ggplot2)
library(ggpubr)

df = read.table("RPKM_static_dynamic.txt",header = TRUE, sep = '\t')

compare_means(df ~ day0,data=df,paired=TRUE)
p <- ggpaired(df,x = "condition",y = "day3",color = "condition", palette = "jco")
p + stat_compare_means(paired = TRUE)

ggplot(data = df, aes(x = condition, y = day0)) +
  geom_tile()

library("GenomicFeatures")
library(GenomicRanges)

chrom.info = read.table(file = "hg19.chrom.sizes_main.txt",header = T)

gencode<-makeTxDbFromGFF("gencode.v26lift37.annotation.gtf",  format="gtf", chrominfo=chrom.info, dataSource=paste("http://genome.ucsc.edu/cgi-bin/hgTables"), organism="Homo sapiens")

gencode_trlengths <- transcriptLengths(gencode,with.cds_len = TRUE,with.utr3_len = TRUE,with.utr5_len = TRUE)

counts<-as.matrix(read.csv("transcript_count_matrix_day_0_3_rep2_rep3.csv",row.names="transcript_id"))

count.df <- read.csv("transcript_count_matrix_day_0_3_rep2_rep3.csv", header = TRUE)

count.df$tx_len = 0
for ( i in 180872:nrow(count.df)){
  x = strsplit(as.character(count.df$transcript_id[i]), split = '_')
  y = strsplit(unlist(x)[1],":")
  if (grepl("ENST",unlist(y)[1], fixed = FALSE)){
    j = grep(pattern = unlist(y)[1], as.character(gencode_trlengths$tx_name))
    if(j!= 0)
      count.df$tx_len[i] = gencode_trlengths$tx_len[j]
  }
  if(grepl("ENST",unlist(y)[2],fixed = FALSE))
    j = grep(pattern = unlist(y)[2], as.character(gencode_trlengths$tx_name))
    if(j!= 0)
      count.df$tx_len[i] = gencode_trlengths$tx_len[j]
  j = 0
  if(i%%10000 == 0)
    print(i)
}

#write.table(count.df,"counts_with_tx_lengths_day0_day3.txt", sep = '\t')



