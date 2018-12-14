#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -J mapBowtie2.job 
#SBATCH -p core -n 6
#SBATCH -t 2-10:00:00
#SBATCH --mail-user pelin.akan@scilifelab.se
#SBATCH --mail-type=ALL

base_file_name=$1
sample_name=$2
echo sample_name $sample_name
echo base_file_name $base_file_name

#SBATCH -o $1.out
#SBATCH -e $1.err

module load bioinfo-tools
module load bowtie2/2.2.9
module load samtools/1.3

#Get current directory
curdir=`readlink -f .`
echo "Current directory is " $curdir
#mkdir $sample_name

export LD_LIBRARY_PATH=/sw/apps/bioinfo/bamtools/2.3.0/milou/lib/:$LD_LIBRARY_PATH
#Pair the Single Mappings
/crex/proj/snic2016-7-108/private/bin/PairSingleEndMappings \
$sample_name/$base_file_name\_1.trunc.bam \
$sample_name/$base_file_name\_2.trunc.bam \
$SNIC_TMP/$base_file_name\.trunc.pair.bam

rm $sample_name/$base_file_name\_1.trunc.bam
rm $sample_name/$base_file_name\_2.trunc.bam

#Filter paired data
/crex/proj/snic2016-7-108/nobackup/External/hicup_v0.5.7/hicup_filter \
--digest /crex/proj/sllstore2017025/private/supplementaryFiles/Digest_hg19_MboI_None_11-39-03_22-01-2016.txt \
--outdir $SNIC_TMP \
$SNIC_TMP/$base_file_name\.trunc.pair.bam

#Remove duplicates
/crex/proj/snic2016-7-108/nobackup/External/hicup_v0.5.7/hicup_deduplicator \
--outdir $SNIC_TMP \
--zip \
$SNIC_TMP/$base_file_name\.trunc.filt.sam

#Sort and index BAM
samtools sort -O bam -o $SNIC_TMP/$base_file_name\.trunc.dedup.sorted.bam -T $SNIC_TMP -@ 2 $SNIC_TMP/$base_file_name\.trunc.dedup.bam
samtools index $SNIC_TMP/$base_file_name\.trunc.dedup.sorted.bam $SNIC_TMP/$base_file_name\.trunc.dedup.sorted.bai

rm $SNIC_TMP/*.sam
rm $SNIC_TMP/*.dedup.bam
rm $SNIC_TMP/*.filt.bam
rm $SNIC_TMP/*.pair.bam
rm $SNIC_TMP/*.trunc.bam
rm $SNIC_TMP/*.fastq.gz
rm $SNIC_TMP/*.fastq
rm $SNIC_TMP/hg19.*
rm -r $SNIC_TMP/hicup_filter_ditag_rejects*
#rm $sample_name/*fastq.gz

echo -ne '\n\nCopying back the resulting analysis files at '
date
rsync -rptoDLv $SNIC_TMP// $curdir/$sample_name//
#rsync -rptoDLv $SNIC_NOBACKUP/$sample_name// $curdir/$sample_name//
#rm -r $SNIC_NOBACKUP/$sample_name//
