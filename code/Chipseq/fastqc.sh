#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 2:00:00
#SBATCH -p core -n 2
#SBATCH -J fastqc
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

module load bioinfo-tools
module load FastQC

mkdir -p /crex/proj/sllstore2017025/nobackup/private/BAV/Chipseq/data/FASTQC

input=/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_54/fastq/
output=/crex/proj/sllstore2017025/nobackup/private/BAV/Chipseq/data/FASTQC

cd ${input}
fastqc -t 2 -f fastq -o ${output} BAV*.gz
