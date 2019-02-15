#!/bin/bash -l

#SBATCH -A g2018023
#SBATCH -o BAV2375_HC_SNP.out
#SBATCH -e BAV2375_HC_SNP.err
#SBATCH -p node -n 16
#SBATCH -t 25:00:00
#SBATCH -J call_snps_BAV2375
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.3
module load picard
module load GATK/3.8-0
module load vcftools
module load vcflib
module load tabix


#add dummy read groups for gatk
java -Xmx8g -jar /sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar AddOrReplaceReadGroups \
I=/crex/proj/sllstore2017025/nobackup/private/BAV/bamfiles/BAV2375/BAV2375.trunc.dedup.sorted.bam \
O=BAV2375.trunc.dedup.sortedbycoord.withreadgroups.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=20 \
VALIDATION_STRINGENCY=SILENT
