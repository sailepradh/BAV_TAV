#Original code by Tobias Everhorn, 2018
## install library before using all of these 

#install.packages("tidyverse")
#install.packages("readr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges", version = "3.8")
#install.packages("tictoc") 


setwd("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/HICCUP_files/")
rm(list=ls())

library(readr)
library(dplyr)
library(GenomicRanges)
library(tictoc)

### Load input file
# Created using the following bedtools line:
# bedtools intersect -wao -a Digest_mm10_MboI_Full.sorted.bed -b All_PeaksMerged_v1_FINAL.txt > Fragments_andOverlapInfo.txt

Frag_overlap_info <- read.delim("Fragments_andOverlapInfo.txt", header = FALSE, stringsAsFactors = FALSE)

#chrM    1       741     1       1       None    Re1     .       -1      -1      0
#chrM    742     952     2       2       Re1     Re1     .       -1      -1      0
#chrM    953     1228    3       3       Re1     Re1     .       -1      -1      0
#chrM    1229    2350    4       4       Re1     Re1     .       -1      -1      0
#chr1    6649876 6650541 9959    9959    Re1     Re1     .       -1      -1      0
#chr1    6650542 6650826 9960    9960    Re1     Re1     .       -1      -1      0
#chr1    6650827 6651467 9961    9961    Re1     Re1     chr1    6651395 6651495 72
#chr1    6651468 6651722 9962    9962    Re1     Re1     chr1    6651395 6651495 27
#chr1    6651723 6652374 9963    9963    Re1     Re1     .       -1      -1      0

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

### Merge fragments within peak regions
new_start=integer(1)
new_end=integer(1)
q=as.integer("0")

library(tictoc)
tic()
i=
j = 1 #Counter for the new digest
while(i < nrow(Frag_overlap_info)) {
  #If fragment has no overlap
  if(Frag_overlap_info[i,11] == 0){
    new_Digest_1[j] <- Frag_overlap_info[i,1]
    new_Digest_2[j] <- Frag_overlap_info[i,2]
    new_Digest_3[j] <- Frag_overlap_info[i,3]
    new_Digest_4[j] <- Frag_overlap_info[i,4]
    new_Digest_5[j] <- Frag_overlap_info[i,5]
    new_Digest_6[j] <- Frag_overlap_info[i,6]
    new_Digest_7[j] <- Frag_overlap_info[i,7]
    i = i + 1
    j = j + 1
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
    new_Digest_1[j] <- Frag_overlap_info[i,1]
    new_Digest_2[j] <- new_start
    new_Digest_3[j] <- new_end
    new_Digest_4[j] <- Frag_overlap_info[i,4]
    new_Digest_5[j] <- Frag_overlap_info[i,5]
    new_Digest_6[j] <- Frag_overlap_info[i,6]
    new_Digest_7[j] <- Frag_overlap_info[i,7]
    j = j + 1
    i = q
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

new_Digest = new_Digest[-(j:i),]

write.table(new_Digest, file="DigestFile_wPeaks_chip-atlas-enhancer-marks_BAVTAV.v2.txt", quote=F, sep="\t", row.names=F, col.names=F)
