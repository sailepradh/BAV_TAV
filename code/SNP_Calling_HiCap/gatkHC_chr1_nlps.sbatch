#!/bin/bash -l

#SBATCH -A g2018023
#SBATCH -o gatkHCchr1nlps.out
#SBATCH -e gatkHCchr1nlps.err
#SBATCH -p node -n 16
#SBATCH -t 20:00:00
#SBATCH -J call_snps_BAV2714
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL


mkdir -p $SNIC_TMP/sail

echo "copy bam files"

rsync -rptoDLv /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/bamfiles/BAV2714EC/BAV2714EC.trunc.dedup.sorted.bam $SNIC_TMP/sail

wait

samplename=$SNIC_TMP/sail/BAV2714EC

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
java -Xmx2g -jar /sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar AddOrReplaceReadGroups.jar \
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
#rsync -rptoDLv $samplename.sorted.dedup.real.recal.bam /pica/w1/snic2016-7-108_nobackup/thp1.wgs/thp1_gatk_snp

wait

echo "call snps"
#Call SNPs
java -Xmx4g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-L chr1 \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta  \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.recal.bam \
-nct 16 \
-pcrModel NONE \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--dbsnp /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf \
-o $samplename.chr1.snp.g.vcf

wait

#copy back snp file
echo "copy back snp files"
rsync -rptoDLv $samplename.chr1.snp.g.vcf /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/BAV2714EC_chr1.vcf
