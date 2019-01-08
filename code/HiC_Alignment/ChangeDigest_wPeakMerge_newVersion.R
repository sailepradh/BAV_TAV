setwd("C:/Users/user/Documents/R/GenerateRandom")
rm(list=ls())
library(readr)
library(dplyr)
library(GenomicRanges)

### Load input file
# Created using the following bedtools line:
# bedtools intersect -wao -a Digest_mm10_MboI_Full.sorted.bed -b All_PeaksMerged_v1_FINAL.txt > Fragments_w_overlap_v8.txt

Frag_overlap_info <- read.delim("Fragments_w_overlap_RANDOM_v8.txt", header = FALSE, stringsAsFactors = FALSE)

#### Define frames and vectors
new_Digest<- data.frame(character(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info)),
                        character(nrow(Frag_overlap_info)),
                        character(nrow(Frag_overlap_info))
)


new_Digest_1=character(nrow(Frag_overlap_info))
new_Digest_2=integer(nrow(Frag_overlap_info))
new_Digest_3=integer(nrow(Frag_overlap_info))
new_Digest_4=integer(nrow(Frag_overlap_info))
new_Digest_5=integer(nrow(Frag_overlap_info))
new_Digest_6=character(nrow(Frag_overlap_info))             
new_Digest_7=character(nrow(Frag_overlap_info))   

merged_Fragments<- data.frame(character(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info)),
                        integer(nrow(Frag_overlap_info))
)

merged_Fragments_1=character(nrow(Frag_overlap_info))
merged_Fragments_2=integer(nrow(Frag_overlap_info))
merged_Fragments_3=integer(nrow(Frag_overlap_info))

### Merge fragments within peak regions
new_start=integer(1)
new_end=integer(1)
q=as.integer("0")

library(tictoc)
tic()
i=1
while(i < nrow(Frag_overlap_info)) {
  #If fragment has no overlap
  if(Frag_overlap_info[i,11] == 0){
    new_Digest_1[i] <- Frag_overlap_info[i,1]
    new_Digest_2[i] <- Frag_overlap_info[i,2]
    new_Digest_3[i] <- Frag_overlap_info[i,3]
    new_Digest_4[i] <- Frag_overlap_info[i,4]
    new_Digest_5[i] <- Frag_overlap_info[i,5]
    new_Digest_6[i] <- Frag_overlap_info[i,6]
    new_Digest_7[i] <- Frag_overlap_info[i,7]
    i=i+1
  #If fragment has overlap  
  }else if(Frag_overlap_info[i,11] != 0 ) {
    new_start=Frag_overlap_info[i,2]
    new_end=Frag_overlap_info[i,3]
    q=i #Index for "connected" fragment
    x=0 #Counter for number of fragments in a streak that are merged
    #Check how many fragments are overlapping to same peak, 10000 =max number of fragments to check forward (large arbitrary)
    #Change largest end position and smallest start position of the fragments that overlap to same peak and merge using the LAST end position
    for(k in 1:10000){
      if(Frag_overlap_info[i+k,11] != 0  & Frag_overlap_info[i,1]==Frag_overlap_info[i+k,1]){
        new_end=Frag_overlap_info[i+k,3]
        q=i+k
        x=x+1
        next()
      }else{
        q=i+k
        break
      }
    }
    new_Digest_1[i] <- Frag_overlap_info[i,1]
    new_Digest_2[i] <- new_start
    new_Digest_3[i] <- new_end
    new_Digest_4[i] <- Frag_overlap_info[i,4]
    new_Digest_5[i] <- Frag_overlap_info[i,5]
    new_Digest_6[i] <- Frag_overlap_info[i,6]
    new_Digest_7[i] <- Frag_overlap_info[i,7]
    
    ## Collect all merged fragments if one wants to look at them
    if(x>0){   
    merged_Fragments_1[i] <- Frag_overlap_info[i,1]
    merged_Fragments_2[i] <- new_start
    merged_Fragments_3[i] <- new_end
    }
    i <- q
  }
}
toc()


# Print While-loop output to one dataframe to create new digest file
new_Digest[,1] = new_Digest_1
new_Digest[,2] = new_Digest_2
new_Digest[,3] = new_Digest_3
new_Digest[,4] = new_Digest_4
new_Digest[,5] = new_Digest_5
new_Digest[,6] = new_Digest_6
new_Digest[,7] = new_Digest_7

names(new_Digest)[1] = "Chromosome"
names(new_Digest)[2] = "Fragment_Start_Position"
names(new_Digest)[3] = "Fragment_End_Position"
names(new_Digest)[4] = "Fragment_Number"
names(new_Digest)[5] = "RE1_Fragment_Number"
names(new_Digest)[6] = "5'_Restriction_Site"
names(new_Digest)[7] = "3'_Restriction_Site"
df = new_Digest

# Print all merged fragments if one wants to look at them
merged_Fragments[,1] = merged_Fragments_1
merged_Fragments[,2] = merged_Fragments_2
merged_Fragments[,3] = merged_Fragments_3
names(merged_Fragments)[1] = "Chromosome"
names(merged_Fragments)[2] = "Fragment_Start_Position"
names(merged_Fragments)[3] = "Fragment_End_Position"

merged_Fragments = merged_Fragments[(merged_Fragments$Fragment_Start_Position > 0),]

###Adjust fragment number , mouse version
chr1=df[(df$Chromosome == "chr1"),]
chr1[,4]=seq(1,nrow(chr1),by=1)
chr1[,5]=seq(1,nrow(chr1),by=1)

chr2=df[(df$Chromosome == "chr2"),]
chr2[,4]=seq(1,nrow(chr2),by=1)
chr2[,5]=seq(1,nrow(chr2),by=1)

chr3=df[(df$Chromosome == "chr3"),]
chr3[,4]=seq(1,nrow(chr3),by=1)
chr3[,5]=seq(1,nrow(chr3),by=1)

chr4=df[(df$Chromosome == "chr4"),]
chr4[,4]=seq(1,nrow(chr4),by=1)
chr4[,5]=seq(1,nrow(chr4),by=1)

chr5=df[(df$Chromosome == "chr5"),]
chr5[,4]=seq(1,nrow(chr5),by=1)
chr5[,5]=seq(1,nrow(chr5),by=1)

chr6=df[(df$Chromosome == "chr6"),]
chr6[,4]=seq(1,nrow(chr6),by=1)
chr6[,5]=seq(1,nrow(chr6),by=1)

chr7=df[(df$Chromosome == "chr7"),]
chr7[,4]=seq(1,nrow(chr7),by=1)
chr7[,5]=seq(1,nrow(chr7),by=1)

chr8=df[(df$Chromosome == "chr8"),]
chr8[,4]=seq(1,nrow(chr8),by=1)
chr8[,5]=seq(1,nrow(chr8),by=1)

chr9=df[(df$Chromosome == "chr9"),]
chr9[,4]=seq(1,nrow(chr9),by=1)
chr9[,5]=seq(1,nrow(chr9),by=1)

chr10=df[(df$Chromosome == "chr10"),]
chr10[,4]=seq(1,nrow(chr10),by=1)
chr10[,5]=seq(1,nrow(chr10),by=1)

chr11=df[(df$Chromosome == "chr11"),]
chr11[,4]=seq(1,nrow(chr11),by=1)
chr11[,5]=seq(1,nrow(chr11),by=1)

chr12=df[(df$Chromosome == "chr12"),]
chr12[,4]=seq(1,nrow(chr12),by=1)
chr12[,5]=seq(1,nrow(chr12),by=1)

chr13=df[(df$Chromosome == "chr13"),]
chr13[,4]=seq(1,nrow(chr13),by=1)
chr13[,5]=seq(1,nrow(chr13),by=1)

chr14=df[(df$Chromosome == "chr14"),]
chr14[,4]=seq(1,nrow(chr14),by=1)
chr14[,5]=seq(1,nrow(chr14),by=1)

chr15=df[(df$Chromosome == "chr15"),]
chr15[,4]=seq(1,nrow(chr15),by=1)
chr15[,5]=seq(1,nrow(chr15),by=1)

chr16=df[(df$Chromosome == "chr16"),]
chr16[,4]=seq(1,nrow(chr16),by=1)
chr16[,5]=seq(1,nrow(chr16),by=1)

chr17=df[(df$Chromosome == "chr17"),]
chr17[,4]=seq(1,nrow(chr17),by=1)
chr17[,5]=seq(1,nrow(chr17),by=1)

chr18=df[(df$Chromosome == "chr18"),]
chr18[,4]=seq(1,nrow(chr18),by=1)
chr18[,5]=seq(1,nrow(chr18),by=1)

chr19=df[(df$Chromosome == "chr19"),]
chr19[,4]=seq(1,nrow(chr19),by=1)
chr19[,5]=seq(1,nrow(chr19),by=1)

chrx=df[(df$Chromosome == "chrX"),]
chrx[,4]=seq(1,nrow(chrx),by=1)
chrx[,5]=seq(1,nrow(chrx),by=1)

chry=df[(df$Chromosome == "chrY"),]
chry[,4]=seq(1,nrow(chry),by=1)
chry[,5]=seq(1,nrow(chry),by=1)

df_final = rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrx,chry)

###Print matrix
library(MASS)
outputMat = df_final
# write.matrix(outputMat,'DigestFile_NEW_wPeaks_adjusted_181115.txt',sep = "\t")
write.matrix(outputMat,'DigestFile_merged_RANDOM_181210.txt',sep = "\t")
write.table(merged_Fragments, file="mergedFragments_RANDOMdigest_181210x.txt", quote=F, sep="\t", row.names=F, col.names=F)
