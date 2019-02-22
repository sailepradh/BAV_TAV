## Generate the random set from Probe -Distal Filtered_dataset
## Written by : Pelin Sahlen
## Recoded by : Sailendra pradhananga
## Email : sailendra.pradhananga@scilifelab.se
## 20/02/2019

## install these libraries for the smooth running
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("HelloRanges", version = "3.8")
# 
# install.packages("bedr")
# install.packages("MASS")
# install.packages("fitdistrplus")

## Loading the library
library(HelloRanges)
library(bedr)
library(MASS)
library(fitdistrplus)

## setting paths for the Random directory
setwd("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_callsFeb/Differential_Interaction/RandomSets/")
Directory = "/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_callsFeb/Differential_Interaction/RandomSets/"

## random sets intialization and directory
set.seed(5)
NofRandomSets = 10
ObservedSet = "sorted_uniq_PD_Interaction"
fname = paste0(Directory,ObservedSet)
fname = paste0(fname,".bed")

## Reading the file
obs = read.table(
                 file = fname,
                 sep = "\t",
                 header = FALSE
                 )

#fit nlog distribution to obtain mu and sigma
x = obs$V3 - obs$V2
fit_ln <- fitdist(x, "lnorm")

#Obtain the genic ratio of observed set
gr_a <- import(fname)
gr_b <- import("hg19.gencode.wholegene.bed")
ans <- subsetByOverlaps(
                        gr_a,
                        gr_b,
                        ignore.strand = TRUE
                        )
genic = length(ans)
intergenic = length(gr_a) - genic

#Generate Random Sets
for (r in 1:NofRandomSets){
  rr1 = get.random.regions(
                           n = nrow(obs)*5,
                           species = "human",
                           build = "hg19",
                           size.mean =  fit_ln$estimate[1],
                           size.sd =  fit_ln$estimate[2]
                          )
  r1 <- rr1[(nchar(rr1$chr)<6),]
  df = data.frame(
                  seqnames = r1$chr,
                  starts = as.integer(r1$start),
                  ends = as.integer(r1$end)
                  )
  write.table(
              df,
              file="foo.bed",
              quote=F,
              sep="\t",
              row.names=F,
              col.names=F
              )

  #Remove those present in observed set
  gr_a <- import("foo.bed")
  gr_b <- import(fname)
  rr2 = subsetByOverlaps(
                         gr_a,
                         gr_b,
                         invert = TRUE,
                         ignore.strand = TRUE
                        )
  df = data.frame(
                  seqnames = seqnames(rr2),
                  starts = start(rr2),
                  ends = end(rr2)
                  )
  write.table(
              df,
              file="foo2.bed",
              quote=F,
              sep="\t",
              row.names=F,
              col.names=F
              )
  #Balance the ratio of intragenic regions
  gr_a <- import("foo2.bed")
  gr_b <- import("hg19.gencode.wholegene.bed")
  rr3 <- subsetByOverlaps(
                          gr_a,
                          gr_b,
                          ignore.strand = TRUE
                          )
  rr3s = rr3[1:genic]
  df = data.frame(
                  seqnames = seqnames(rr3s),
                  starts = start(rr3s),
                  ends = end(rr3s)
                  )
  fn = paste(ObservedSet,"_",as.character(r),".bed")
  fn = gsub(" ", "",fn)
  write.table(
              df,
              file=fn,
              quote=F,
              sep="\t",
              row.names=F, 
              col.names=F
              )
  rr4 <- subsetByOverlaps(
                          gr_a,
                          gr_b,invert = TRUE,
                          ignore.strand = TRUE
                          )
  rr4s <- rr4[1:intergenic]
  df = data.frame(
                  seqnames = seqnames(rr4s),
                  starts = start(rr4s),
                  ends = end(rr4s)
                  )
  write.table(
              df,
              file=fn,
              quote=F,
              sep="\t",
              row.names=F,
              col.names=F,
              append = TRUE
              )
  print(r)
}
