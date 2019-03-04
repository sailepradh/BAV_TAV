#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -J mapBowtie02_02_27.job
#SBATCH -p node
#SBATCH -t 1-12:00:00
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
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

#Create Symbolic Link to genome
ln -s /sw/data/uppnex/reference/biodata/genomes/Hsapiens/hg19/bowtie2/hg19* $SNIC_TMP/.

echo "Map with Bowtie 2" *.trunc.fastq.gz
#time bowtie2 -p 10 --very-sensitive --reorder -x $SNIC_TMP/hg19 -U $SNIC_NOBACKUP/$sample_name/$base_file_name\_R1_001.trunc.fastq.gz >$SNIC_TMP/$base_file_name\_1.trunc.sam &
#time bowtie2 -p 10 --very-sensitive --reorder -x $SNIC_TMP/hg19 -U $SNIC_NOBACKUP/$sample_name/$base_file_name\_R2_001.trunc.fastq.gz >$SNIC_TMP/$base_file_name\_2.trunc.sam &

time bowtie2 \
-p 10 \
--very-sensitive \
--reorder \
-x $SNIC_TMP/hg19 \
-U $sample_name/$base_file_name\_R1_001.trunc.fastq.gz \
| samtools view -bS - > $sample_name/$base_file_name\_1.trunc.bam &


time bowtie2 \
-p 10 \
--very-sensitive \
--reorder \
-x $SNIC_TMP/hg19 \
-U $sample_name/$base_file_name\_R2_001.trunc.fastq.gz \
| samtools view -bS - > $sample_name/$base_file_name\_2.trunc.bam &

wait

rm $sample_name/*trunc.fastq.gz

#samtools view -bS $SNIC_TMP/$base_file_name\_1.trunc.sam >$SNIC_NOBACKUP/$sample_name/$base_file_name\_1.trunc.bam &
#samtools view -bS $SNIC_TMP/$base_file_name\_2.trunc.sam >$SNIC_NOBACKUP/$sample_name/$base_file_name\_2.trunc.bam &

#rsync -rptoDLv $SNIC_NOBACKUP/$sample_name/$base_file_name\_1.trunc.bam /proj/sllstore2017025/nobackup/private/TAV/bamfiles/updated/ &
#rsync -rptoDLv $SNIC_NOBACKUP/$sample_name/$base_file_name\_2.trunc.bam /proj/sllstore2017025/nobackup/private/TAV/bamfiles/updated/ &

wait

sbatch 03_Process_BAM.sbatch $base_file_name $sample_name
