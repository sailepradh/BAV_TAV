#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -J mapBowtie2.job 
#SBATCH -p core -n 6
#SBATCH -t 1-12:00:00
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

mkdir -p $sample_name
echo "Created a directory within current directory " $sample_name

#Merge FASTQ files
cat \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_14/BAV2375rep1de_S1_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_14/BAV2375rep2de_S2_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_15/BAV2375rep1de_S2_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_15/BAV2375rep2de_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_16/BAV2375rep1de_S2_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_17/BAV2375rep1des4_S5_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_17/BAV2375rep2des4_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_22/BAV2375rep2des4_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_23/BAV2375rep1des4_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_23/BAV2375rep2des4_S2_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_24/BAV2375rep1des4_S4_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_24/BAV2375rep2des4_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_25/BAV2375rep1des4_S6_R1_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_25/BAV2375rep2des4_S5_R1_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R1_001.fastq.gz &

cat \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_14/BAV2375rep1de_S1_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_14/BAV2375rep2de_S2_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_15/BAV2375rep1de_S2_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_15/BAV2375rep2de_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_16/BAV2375rep1de_S2_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_17/BAV2375rep1des4_S5_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_17/BAV2375rep2des4_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_22/BAV2375rep2des4_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_23/BAV2375rep1des4_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_23/BAV2375rep2des4_S2_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_24/BAV2375rep1des4_S4_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_24/BAV2375rep2des4_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_25/BAV2375rep1des4_S6_R2_001.fastq.gz \
/crex/proj/sllstore2017025/nobackup/private/BAV/BAV_fastq_files/NextSeq_25/BAV2375rep2des4_S5_R2_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R2_001.fastq.gz &

wait

#Hicup_truncator
/crex/proj/sllstore2017025/nobackup/private/External/hicup_v0.5.7/hicup_truncater \
--zip \
--threads 2 \
--re1 ^GATC \
--outdir $sample_name/ \
$SNIC_TMP/$base_file_name\_R1_001.fastq.gz $SNIC_TMP/$base_file_name\_R2_001.fastq.gz

sbatch 02_MapPairs.sbatch $base_file_name $sample_name
