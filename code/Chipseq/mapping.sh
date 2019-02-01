#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 3:00:00
#SBATCH -p core -n 10
#SBATCH -J mapping$1
#SBATCH -o mapping.out
#SBATCH -e mapping.err
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

module load bioinfo-tools
module load bowtie2
module load samtools

input=/crex/proj/sllstore2017025/private/SEQUENCING_RUNS_FASTQ/NextSeq_54/fastq/
cd ${input}

base_file_name=$1
name=$base_file_name"_R1_001.fastq.gz"


bowtie2 \
    -x /sw/data/uppnex/reference/Homo_sapiens/GRCh37/program_files/bowtie2/concat \
    -U $name \
    -p 10 1> $base_file_name".sam" 2> $base_file_name".bt2.log"

wait

samtools view -bhS $base_file_name".sam" > $base_file_name".bam"
