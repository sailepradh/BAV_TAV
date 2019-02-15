#!/bin/bash -l

#SBATCH -A a2014205
#SBATCH -o gatkHCchr1wlps.out
#SBATCH -e gatkHCchr1wlps.err
#SBATCH -p node -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J call_snpswLPSchr1
#SBATCH --mail-user aanil@scilifelab.se
#SBATCH --mail-type=ALL


mkdir -p $SNIC_TMP/anan/

echo "copy bam files"

rsync -rptoDLv /proj/b2015060/nobackup/ANALYSIS/THP1_merged/THP1.wLPS.rep$1/THP1.wLPS.rep$1.trunc.dedup.sortedbycoord.bam $SNIC_TMP/anan

wait

samplename=$SNIC_TMP/anan/THP1.wLPS.rep$1


rsync -rptoDLv /pica/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf $SNIC_TMP/anan


module load bioinfo-tools
module load samtools/1.3
module load picard
module load GATK/3.2.2
module load vcftools
module load vcflib
module load tabix

bgzip $SNIC_TMP/anan/dbsnp_138.hg19.vcf
tabix -p vcf $SNIC_TMP/anan/dbsnp_138.hg19.vcf.gz


#add dummy read groups for gatk
java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar \
I=$samplename.trunc.dedup.sortedbycoord.bam \
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
java -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /pica/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.bam \
-o $samplename.recal_data.table \
-cov 'ReadGroupCovariate' \
-cov 'QualityScoreCovariate' \
-cov 'CycleCovariate' \
-cov 'ContextCovariate' \
-knownSites $SNIC_TMP/anan/dbsnp_138.hg19.vcf.gz \
-nct 16

echo "copy back quality score table"
rsync -rptoDLv $samplename.recal_data.table /proj/b2015060/nobackup/ANALYSIS/THP1_merged/SNPvcf

wait

#Print Reads
echo "print reads"
java -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
-T PrintReads \
-R /pica/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
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
java -Xmx4g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-L chr1 \
-R /pica/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-I $samplename.trunc.dedup.sortedbycoord.withreadgroups.recal.bam \
-nct 16 \
-pcrModel NONE \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--dbsnp /pica/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf \
-o $samplename.chr1.snp.g.vcf

wait

#copy back snp file
echo "copy back snp files"
rsync -rptoDLv $samplename.chr1.snp.g.vcf /proj/b2015060/nobackup/ANALYSIS/THP1_merged/SNPvcf
