#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
An updated script on variant manipulator from the HiCap data.
Rare variant are 0.5% within a population at the current criteria.

current use:
python3 ./HiCapVariantparser.py
/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected/test_2500.txt
/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/SNP_calling/Individual_and_combined_VCFs/April_05/Annotated_snp_indel.vcf.gz
--o /Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected/SNP_Interaction.txt
'''


import sys
import argparse
import re
from pysam import VariantFile

pattern = re.compile('SwedFreq.AF')


def int2binary(interaction_sample):
    interaction_binary = []
    for int_sta in interaction_sample:
        if int(int_sta[0]) >= 4 and float(int_sta[1]) < 0.01:
            interaction_binary.append("1")
        else:
            interaction_binary.append("0")
    return (interaction_binary)


def Main ():
    parser = argparse.ArgumentParser(description="loading vcf and interaction files")
    parser.add_argument("interactionfile", help = "Interaction calls from HiCap method")
    parser.add_argument("vcfile", help = "Variant calls from either HiCap or sequencing samples")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()
    Vcfin = VariantFile(args.vcfile)

    result_title = ["RefSeqName","TranscriptName","Feature_ID",
                    "Feature_Chr","Feature_Start","Feature_End",
                    "Annotation","Strand","Interactor_Chr","Interactor_Start","Interactor_End","Distance",
                    "SNPs","SNP_ID","Ind_count","Swed_Freq",
                    "TAV2431","TAV2515","TAV2709","BAV2375","BAV2424","BAV2714"]

    with open (args.output, "w") as output_file:
        output_file.write ("\t".join(result_title)+"\n")

    with open (args.interactionfile, 'r') as f:
        next(f)

        for line in f:
            line = line.strip().split("\t")
            all_fields = line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11]
            chr = ((line[8])[3:],line[9], line[10])

            TAV2431 = [line [12],line[13]]
            TAV2515 = [line [15],line[16]]
            TAV2709 = [line [18],line[19]]
            BAV2375 = [line [21],line[22]]
            BAV2424 = [line [24],line[25]]
            BAV2714 = [line [27],line[28]]

            interaction_sample = [TAV2431, TAV2515, TAV2709, BAV2375, BAV2424, BAV2714]
            interaction_binary = int2binary(interaction_sample)

            sample_list = [3,4,5,0,1,2]
            for rec in Vcfin.fetch (chr[0], int(chr[1]), int(chr[2])):
                genotype_binary = []
                for test in rec.samples.values():

                    genotype =  "/".join([str(x) for x in test["GT"]])
                    if genotype == "None/None":
                        continue
                    elif genotype == "0/1" or genotype == "1/1":
                        genotype_binary.append("1")
                    elif genotype == "0/0":
                        genotype_binary.append("0")

                    swed_freq = "0"
                    for f,v in rec.info.iteritems():
                        if pattern.match (f):
                            swed_freq = v

                    if rec.id == None:
                        rec.id = "X"

                sorted_genotype = [x for _,x in sorted(zip(sample_list,genotype_binary))]
                zip_array = list(zip (interaction_binary, sorted_genotype))

                count_0 = 0
                count_1 = 0
                mismatch_int = 0
                mismatch_var = 0

                for a,b in zip_array:
                    if a == '0' and b == '1':
                        mismatch_int = mismatch_int + 1

                    if a == '1'  and b == '0':
                        mismatch_var = mismatch_var + 1

                    if a == '0' and b == '0':
                        count_0 =count_0 + 1

                    if a == '1' and b == '1':
                        count_1 = count_1 + 1

                count = count_0+count_1

                if count == 6:
                    allele = "|".join(rec.alleles)
                    count_int_allele = 0
                    for a,b in zip_array:
                        if (a,b) == ('1','1'):
                            count_int_allele = count_int_allele + 1

                    changed_freq = "".join(str(x) for x in swed_freq)
                    unzip_array =  ["|".join(x) for x in zip_array]
                    snp = (line[8],rec.start,rec.stop,allele,rec.filter.keys()[0])
                    str_snp = "_".join(str(x) for x in snp)

                    result = "\t".join(all_fields),str_snp,rec.id,count_int_allele,changed_freq,"\t".join(unzip_array)
                    combined_result = "\t".join(str(x) for x in result)

                    with open (args.output, "a") as output_file:
                        output_file.write (combined_result+"\n")

                if mismatch_int == 1 and mismatch_var == 1 and count_0 == 4 :
                    allele = "|".join(rec.alleles)
                    count_int_allele = 0
                    for a,b in zip_array:
                        if (a,b) == ('0','1'):
                            count_int_allele = count_int_allele + 1


                    changed_freq = "".join(str(x) for x in swed_freq)
                    unzip_array =  ["|".join(x) for x in zip_array]
                    snp = (line[8],rec.start,rec.stop,allele,rec.filter.keys()[0])
                    str_snp = "_".join(str(x) for x in snp)

                    result = "\t".join(all_fields),str_snp,rec.id,count_int_allele,changed_freq,"\t".join(unzip_array)
                    combined_result = "\t".join(str(x) for x in result)

                    with open (args.output, "a") as output_file:
                        output_file.write (combined_result+"\n")

if __name__=="__main__":
    Main()
