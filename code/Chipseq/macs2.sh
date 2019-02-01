#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 1:00:00
#SBATCH -p core -n 1
#SBATCH -J macs2_callpeak
#SBATCH -o macs2_callpeak.out
#SBATCH -e macs2_callpeak.err
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=FAIL
#SBATCH --tmp=20480

module load bioinfo-tools
module load MACS/2.1.0


input=/crex/proj/sllstore2017025/nobackup/private/BAV/Chipseq/data/BAM/
output=/crex/proj/sllstore2017025/nobackup/private/BAV/Chipseq/data/MACS

cd ${input}
echo "Calling peaks for " $1
echo "Control file " $2
macs2 callpeak -t $1 -c $2 -n $1 -f BAM -g hs -q 0.01 --outdir ${output}
echo "Done"
