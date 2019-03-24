#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -J mapBowtie01_BAV2742_03_20.job
#SBATCH -p core -n 6
#SBATCH -t 1-12:00:00
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL

base_file_name=$1
#BAV2764smc_S2
#TAV2675smc_S1
#TAV2695smc_S5

sample_name=$2
#BAV2764smc
#TAV2675smc
#TAV2695smc

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
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_62/fastq/BAV2742SMCdbenr_S1_R1_001.fastq.gz \
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_65/fastq/BAV2742smc_S3_R1_001.fastq.gz \
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_66/fastq/BAV2742smoothm_S1_R1_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R1_001.fastq.gz &

cat \
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_62/fastq/BAV2742SMCdbenr_S1_R2_001.fastq.gz  \
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_65/fastq/BAV2742smc_S3_R2_001.fastq.gz \
/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_66/fastq/BAV2742smoothm_S1_R2_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R2_001.fastq.gz &

#copy the fastq files to the scratch

#cp /crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_65/fastq/$base_file_name\_R1_001.fastq.gz $SNIC_TMP/$base_file_name\_R1_001.fastq.gz &
#cp /crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_65/fastq/$base_file_name\_R2_001.fastq.gz $SNIC_TMP/$base_file_name\_R2_001.fastq.gz &

wait

#Hicup_truncator
/crex/proj/sllstore2017025/nobackup/private/External/hicup_v0.5.7/hicup_truncater \
--zip \
--threads 2 \
--re1 ^GATC \
--outdir $sample_name/ \
$SNIC_TMP/$base_file_name\_R1_001.fastq.gz $SNIC_TMP/$base_file_name\_R2_001.fastq.gz

sbatch 02_MapPairs.sbatch $base_file_name $sample_name
