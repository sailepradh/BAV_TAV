#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p devcore -n 4
#SBATCH -t 1:00:00
#SBATCH -J filterPD
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480


#awk -F '\t' '{if (NR ==1) print $0 ; if ($12 != -1 && $13 >= 4 && $14 < 0.001 && $16 >= 4 && $17 < 0.001 ) print $0}' TAVrun.hg19.Proximities.Probe_Distal.txt >  TAVrun.hg19.Proximities.Probe_Distal_SP4_p0.001_filtered.txt

#awk -F '\t' '{ if (NR ==1) print $0 ; if ($17 != -1 && $18 >= 4 && $19 < 0.001 && $21 >= 4 && $22 < 0.001 ) print $0}' TAVrun.hg19.Proximities.Probe_Probe.txt > TAVrun.hg19.Proximities.Probe_Probe_SP4_p0.001_filtered.txt

#awk -F '\t' '{ if (NR ==1) print $0 ; if ($12 != -1 && $13 >= 4 && $14 < 0.05 && $16 >= 4 && $17 < 0.05 ) print $0}' TAVrun.hg19.Proximities.Probe_Distal.txt > TAVrun.hg19.Proximities.Probe_Distal_SP4_p0.05_filtered.txt


#awk -F '\t' '{ if (NR ==1) print $0 ; if ($17 != -1 && $18 >= 4 && $19 < 0.05 && $21 >= 4 && $22 < 0.05 ) print $0}' TAVrun.hg19.Proximities.Probe_Probe.txt > TAVrun.hg19.Proximities.Probe_Probe_SP4_p0.05_filtered.txt


## Jan-09-2019

awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 1) || ($12 != -1 && $16 >= 1)) print $0}' ../BAV_TAV.Proximities.Probe_Distal.txt  > BAV_TAV.Proximities.Probe_Distal_SP1_filtered_in_one.txt

######## Strict filter

#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.001) || ($12 != -1 && $16 >= 4 && $17 < 0.001)) print $0}'   BAV_TAV.Proximities.Probe_Distal.txt > BAV_TAV.Proximities.Probe_Distal_SP4_p001_filtered.txt


#awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.001) || ($17 != -1 && $21 >= 4 && $22 < 0.001)) print $0}' BAV_TAV.Proximities.Probe_Probe.txt > BAV_TAV.Proximities.Probe_probe_SP4_p001_filtered.txt


######## Lenient filter

#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.05) || ($12 != -1 && $16 >= 4 && $17 < 0.05)) print $0}'   BAV_TAV.Proximities.Probe_Distal.txt > BAV_TAV.Proximities.Probe_Distal_SP4_p05_filtered.txt


#awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.05) || ($17 != -1 && $21 >= 4 && $22 < 0.05)) print $0}' BAV_TAV.Proximities.Probe_Probe.txt > BAV_TAV.Proximities.Probe_probe_SP4_p05_filtered.txt





