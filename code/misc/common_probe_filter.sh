#!/bin/bash -l
#SBATCH -A snic2017-7-367
#SBATCH -p core -n 4
#SBATCH -t 5:00:00
#SBATCH -J mergefiles
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

cd /crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/run2/data/script

python3 test.py
