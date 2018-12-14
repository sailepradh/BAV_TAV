#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -J mapBowtie2.job 
#SBATCH -p core -n 4
#SBATCH -t 1-00:00:00
#SBATCH --mail-user aanil@scilifelab.se
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

mkdir -p $SNIC_NOBACKUP/$sample_name

#Merge FASTQ files
cat \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_34/fastq/TAV2515rep1_S2_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_34/fastq/TAV2515rep2_S4_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_35/fastq/TAV2515rep1_S2_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_35/fastq/TAV2515rep2_S4_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_37/fastq/TAV2515rep1_S2_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_37/fastq/TAV2515rep2_S4_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_42/fastq/TAV2515r1_S3_R1_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_42/fastq/TAV2515r2_S1_R1_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R1_001.fastq.gz &

cat \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_34/fastq/TAV2515rep1_S2_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_34/fastq/TAV2515rep2_S4_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_35/fastq/TAV2515rep1_S2_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_35/fastq/TAV2515rep2_S4_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_37/fastq/TAV2515rep1_S2_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_37/fastq/TAV2515rep2_S4_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_42/fastq/TAV2515r1_S3_R2_001.fastq.gz \
/crex2/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_42/fastq/TAV2515r2_S1_R2_001.fastq.gz \
>$SNIC_TMP/$base_file_name\_R2_001.fastq.gz &

wait
#Hicup_truncater 
/proj/snic2016-7-108/nobackup/External/hicup_v0.5.7/hicup_truncater \
--zip \
--threads 2 \
--re1 ^GATC \
--outdir $SNIC_NOBACKUP/$sample_name/ \
$SNIC_TMP/$base_file_name\_R1_001.fastq.gz $SNIC_TMP/$base_file_name\_R2_001.fastq.gz

sbatch 02_MapPairs.sbatch $base_file_name $sample_name
