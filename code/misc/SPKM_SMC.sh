#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J SPKM_TAV
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480


#awk -F '\t' '{ if (($12 != -1 && $13 >= 4 && $14 < 0.001 && $16 >= 4 && $17 < 0.001) || ($12 != -1 && $19 >= 4 && $20 < 0.001 && $22 >= 4 && $23 < 0.001) ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt >  THP1run.hg19.Proximities.Probe_Distal_SP4_p0.001_filtered.txt


#awk -F '\t' '{ if (($17 != -1 && $18 >= 4 && $19 < 0.001 && $21 >= 4 && $22 < 0.001) || ($17 != -1 && $24 >= 4 && $25 < 0.001 && $27 >= 4 && $28 < 0.001) ) print $0}' THP1run.hg19.Proximities.Probe_Probe.txt > THP1run.hg19.Proximities.Probe_Probe_SP4_p0.001_filtered.txt

#awk -F '\t' '{ if (($12 != -1 && $13 >= 4 && $14 < 0.05 && $16 >= 4 && $17 < 0.05) || ($12 != -1 && $19 >= 4 && $20 < 0.05 && $22 >= 4 && $23 < 0.05) ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt > THP1run.hg19.Proximities.Probe_Distal_SP4_p0.05_filtered.txt


#awk -F '\t' '{ if (($17 != -1 && $18 >= 4 && $19 < 0.05 && $21 >= 4 && $22 < 0.05) || ($17 != -1 && $24 >= 4 && $25 < 0.05 && $27 >= 4 && $28 < 0.05) ) print $0}' THP1run.hg19.Proximities.Probe_Probe.txt > THP1run.hg19.Proximities.Probe_Probe_SP4_p0.05_filtered.txt

#awk -F '\t' '{ if ($12 != -1 ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt >  THP1run.hg19.Proximities.Probe_Distal_intrachrom_filtered.txt

#wc -l THP1run.hg19.Proximities.Probe_Distal_intrachrom_filtered.txt

#awk -F '\t' -v OFS="\t" '{ if (($12 != -1 && $13 >= 1) || ($12 != -1 && $16 >= 1))  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' THP1run.hg19.Proximities.Probe_Distal.txt >  tmp1

#wc -l tmp1

#awk -F '\t' -v OFS="\t" '{ if (($12 != -1 && $19 >= 1) || ($12 != -1 && $22 >= 1))  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$19,$20,$21,$22,$23,$24}' THP1run.hg19.Proximities.Probe_Distal.txt > tmp2

#wc -l tmp2

#awk -F '\t' -v OFS="\t" '{ if (($17 != -1 && $18 >= 1) || ($17 != -1 && $21 >= 1)) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14,$17,$18,$19,$20,$21,$22,$23}' THP1run.hg19.Proximities.Probe_Probe.txt | sed 1d > tmp3

#wc -l tmp3

#awk -F '\t' -v OFS="\t" '{ if (($17 != -1 && $24 >= 1) || ($17 != -1 && $27 >= 1)) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14,$17,$24,$25,$26,$27,$28,$29}' THP1run.hg19.Proximities.Probe_Probe.txt | sed 1d  > tmp4

#wc -l tmp4

#echo "temp files created"

#cat tmp1 tmp3 > THP1run.hg19.Proximities.both_nLPS.txt
#rm -r tmp1 tmp3
#cat tmp2 tmp4 > THP1run.hg19.Proximities.both_wLPS.txt
#rm -r tmp2 tmp4

#awk '{print $1}' THP1run.hg19.Proximities.both_nLPS.txt | sort | uniq > features_nLPS.txt
#awk '{print $1}' THP1run.hg19.Proximities.both_wLPS.txt | sort | uniq > features_wLPS.txt

#echo "all done"

#time python3 hicap_expr_PP_raw_18_Sep_2018.py SMC_filtered_16_11_18.txt  Common_Expressed_gene_18_Nov.txt  --o SMC_CPM_Exprs_18_Nov.txt
#time python3 hicap_expr_PP_raw_18_Sep_2018.py SMC_filtered_16_11_18.txt  Common_UnExpressed_gene_18_Nov.txt --o SMC_CPM_Unexprs_18_Nov.txt

#hicap_expr_PP_raw_SMC_21_Nov.py
time python3 hicap_expr_PP_raw_TAV.py TAV_All_Filtered_three_replicates.txt Genes_Target_BAV_TAV.txt  --o TAV_all_SP.txt
