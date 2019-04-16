## Rare and common vairant interpretation from the interaction data 
## This is the down-streaming analysis of the variants and interaction table generated from the earlier python script. The aim here is to visualize the result and process for further analysis.
## The question I am gt

library(reshape2)

setwd ("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected")

SNP_interaction_table = read.table("SNP_Interaction.txt",
                                   sep="\t", header=TRUE,
                                   stringsAsFactors=F)
head (SNP_interaction_table)
gene_enhancer <- as.data.frame (table (SNP_interaction_table$RefSeqName))

head (gene_enhancer)
sorted_table = gene_enhancer[order(-gene_enhancer$Freq), ]
head (sorted_table)

lst_sample = split(SNP_interaction_table, SNP_interaction_table$Ind_count)
 
rare = function(x) {
  y <- subset(x, x$Swed_Freq < 0.005)
  dim_y = dim(y)
  return(dim_y)
}

lowfreq = function(x) {
  y <- subset(x, x$Swed_Freq <= 0.01 &  x$Swed_Freq >= 0.005 )
  dim_y = dim(y)
  return(dim_y)
}

common = function(x) {
  y <- subset(x, x$Swed_Freq > 0.01)
  dim_y = dim(y)
  return(dim_y)
}

lst_rare <- lapply (lst_sample, function(x) rare(x))
rare_num = do.call (rbind, lst_rare)
                    
lst_lf <- lapply(lst_sample, function(x) lowfreq(x))
lowfreq_num = do.call (rbind,lst_lf)

lst_comm <- lapply (lst_sample, function(x) common(x))
common_num = do.call (rbind,lst_comm)

combined_data = cbind(rare_num[,1],
lowfreq_num[,1],
common_num[,1])

rownames(combined_data)
colnames(combined_data) <-c("Rare", "LowFreq", "Common")
combined_df <- as.data.frame(combined_data)

combined_df$category <- row.names(combined_df)
mdfr <- melt(combined_df , id.vars = "category")

library(scales)
library(ggplot2)

bp<- ggplot(mdfr, aes(x="", value, fill = variable))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie

(p <- ggplot(mdfr, aes(category, value, fill = variable)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = number)
)


