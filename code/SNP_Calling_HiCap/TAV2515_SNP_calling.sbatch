#!/bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -o TAV2515_HC_SNP.out
#SBATCH -e TAV2515_HC_SNP.err
#SBATCH -p node -n 16
#SBATCH -t 22:00:00
#SBATCH -J call_snps_TAV2515
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL

mkdir -p $SNIC_TMP/sail

echo "copy bam files"

rsync -rptoDLv /crex/proj/sllstore2017025/nobackup/private/TAV/bamfiles/updated/TAV2515/TAV2515.trunc.dedup.sorted.bam $SNIC_TMP/sail

wait

samplename=$SNIC_TMP/sail/TAV2515

rsync -rptoDLv /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf $SNIC_TMP/sail

module load bioinfo-tools
module load samtools/1.3
module load picard
module load GATK/3.8-0
module load vcftools
module load vcflib
module load tabix

bgzip $SNIC_TMP/sail/dbsnp_138.hg19.vcf
tabix -p vcf $SNIC_TMP/sail/dbsnp_138.hg19.vcf.gz

wait

#add dummy read groups for gatk
java -Xmx4g -jar /sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar AddOrReplaceReadGroups \
I=$samplename.trunc.dedup.sorted.bam \
O=$samplename.trunc.dedup.sortedbycoord.withreadgroups.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=20 \
VALIDATION_STRINGENCY=SILENT

wait

echo "generate bam index files"
samtools index -b $samplename.trunc.dedup.sortedbycoord.withreadgroups.bam

#ReCalibrate Quality Scores
echo "recalibrate quality scores"
java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.bam \
-o $samplename.recal_data.table \
-cov 'ReadGroupCovariate' \
-cov 'QualityScoreCovariate' \
-cov 'CycleCovariate' \
-cov 'ContextCovariate' \
-knownSites $SNIC_TMP/sail/dbsnp_138.hg19.vcf.gz \
-nct 16

echo "copy back quality score table"
rsync -rptoDLv $samplename.recal_data.table /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/

wait

#Print Reads
echo "print reads"
java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.bam \
-baq 'CALCULATE_AS_NECESSARY' \
-BQSR $samplename.recal_data.table \
-o $samplename.trunc.dedup.sortedbycoord.withreadgroups.recal.bam \
-nct 16

echo "copy back recalibrated bam file"
rsync -rptoDLv $samplename.trunc.dedup.sortedbycoord.withreadgroups.recal.bam  /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/TAV2515.trunc.dedup.sortedbycoord.withreadgroups.recal.bam
wait

echo "call snps"
#Call SNPs
java -Xmx4g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta  \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.recal.bam \
-ERC GVCF \
-G Standard \
-G AS_Standard \
-nct 16 \
--dbsnp /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf \
-o $samplename.snp.g.vcf

wait

#copy back snp file
echo "copy back snp files"
rsync -rptoDLv $samplename.snp.g.vcf /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/TAV2515.vcf


echo "VCF file made"
