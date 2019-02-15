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


######## Strict filter
# awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.001 && $16 >= 4 && $17 < 0.001 && $19 >= 4 && $20 < 0.001) || ($12 != -1 && $22 >= 4 && $23 < 0.001 && $25 >= 4 && $26 < 0.001 && $28 >= 4 && $29 < 0.001 )) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP4_p001_filtered.txt

# awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.001 && $21 >= 4 && $22 < 0.001 && $24 >= 4 && $25 < 0.001) || ($17 != -1 && $27 >= 4 && $28 < 0.001 && $30 >= 4 && $31 < 0.001 && $33 >= 4 && $34 < 0.001)) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Probe.txt  > BAVTAV.Proximities.Probe_probe_SP4_p001_filtered.txt


####### lenient filter

#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.01 && $16 >= 4 && $17 < 0.01 && $19 >= 4 && $20 < 0.01) || ($12 != -1 && $22 >= 4 && $23 < 0.01 && $25 >= 4 && $26 < 0.01 && $28 >= 4 && $29 < 0.01 )) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP4_p01_filtered.txt

awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.01 && $21 >= 4 && $22 < 0.01 && $24 >= 4 && $25 < 0.01) || ($17 != -1 && $27 >= 4 && $28 < 0.01 && $30 >= 4 && $31 < 0.01 && $33 >= 4 && $34 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Probe.txt  > BAVTAV.Proximities.Probe_probe_SP4_p01_filtered.txt


##### medium filter

#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.1 && $16 >= 4 && $17 < 0.1 && $19 >= 4 && $20 < 0.1) || ($12 != -1 && $22 >= 4 && $23 < 0.1 && $25 >= 4 && $26 < 0.1 && $28 >= 4 && $29 < 0.1 )) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP4_p1_filtered.txt

#awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.1 && $21 >= 4 && $22 < 0.1 && $24 >= 4 && $25 < 0.1) || ($17 != -1 && $27 >= 4 && $28 < 0.1 && $30 >= 4 && $31 < 0.1 && $33 >= 4 && $34 < 0.1)) print $0}' BAVTAV_mergedDigestPeak.Proximities.Probe_Probe.txt  > BAVTAV.Proximities.Probe_probe_SP4_p1_filtered.txt










## Jan-09-2019

## awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 1) || ($12 != -1 && $16 >= 1) || ($12 != -1 && $19>= 1)  || ( $12 != -1 && $22 >= 1)) print $0}' ../BAVTAV.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP1_filtered_in_one.txt

## individual_filter

#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.001) || ($12 != -1 && $16 >= 4 && $17 < 0.001) || ($12 != -1 && $19>= 4 && $20 < 0.001)  || ( $12 != -1 && $22 >= 4 && $23 < 0.001)) print $0}' ../BAVTAV.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP4_p001_filtered_in_one.txt


#awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.001) || ($17 != -1 && $21 >= 4 && $22 < 0.001) || ($17 != -1 && $24 >= 4 && $25 < 0.001 ) || ($17 != -1 && $27 >= 4 && $28 < 0.001)) print $0}' ../BAVTAV.Proximities.Probe_Probe.txt  > BAVTAV.Proximities.Probe_probe_SP4_p001_filtered_in_one.txt



######## Lenient filter
#awk -F '\t' '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.05 && $16 >= 4 && $17 < 0.05) || ($12 != -1 && $19 >= 4 && $20 < 0.05 && $22 >= 4 && $23 < 0.05)) print $0}' BAVTAV.Proximities.Probe_Distal.txt  > BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered.txt

#awk -F '\t' '{ if (NR ==1) print $0 ; if (($17 != -1 && $18 >= 4 && $19 < 0.05 && $21 >= 4 && $22 < 0.05) || ($17 != -1 && $24 >= 4 && $25 < 0.05 && $27 >= 4 && $28 < 0.05)) print $0}' BAVTAV.Proximities.Probe_Probe.txt  > BAVTAV.Proximities.Probe_probe_SP4_p05_filtered.txt





