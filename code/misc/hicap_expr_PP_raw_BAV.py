#!/usr/bin/env python
# -*- coding: utf-8 -*


"""
python hicap_expr_PP_raw.py /Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/PP/TAVrun.hg19.Proximities.Probe_probe_raw_filtered.txt /Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/Common_Expr_probe_gene.txt -o ../data/raw_external/TAV_PP_Exp_SP_FPKM
"""
import os
import sys
import argparse

def expr_gene(Common_expr_probe):
    genes =[]
    for lines in Common_expr_probe:
        line = lines.strip()
        genes.append(line)
    return(genes)

def int_dict_SP(int_files,expr_genes):

    """
        This function takes the interaction file of PD interaction and returns dictionary of gene name as
        the key and transcript_Id, distal chromosome position, as the list.

        Also it gives the dictionary of supporting pairs

        file = interaction file from HiCap runs with RefSeqName as the gene_ID, TranscriptName as transcript_ID

    """

    results_interaction_SP = {}
    counts_rep1 = 0
    counts_rep2 = 0
    counts_rep3 = 0

    for lines in int_files:
        line = lines.strip()
        fields = line.split("\t")
        counts_rep1 = int(fields[12])+ int(counts_rep1)
        counts_rep2 = int(fields[15])+ int(counts_rep2)
        counts_rep3 = int(fields[18])+ int(counts_rep2)


        if fields[0] in expr_genes:
            interaction_genes = fields[0]
            length = abs(int(fields[9])- int(fields[10]))
            interaction_status_SP = results_interaction_SP.get(interaction_genes,[])
            interaction_status_SP.append( [fields[12], fields[15],field[18],length])
            results_interaction_SP[interaction_genes] = interaction_status_SP
        else:
            pass

    print ("PD_dictionary_done")
    print (counts_rep1,counts_rep2,counts_rep3)
    return (results_interaction_SP, counts_rep1,counts_rep2,counts_rep3)

def SP_CPM(results_interaction, supp_count_rep1,supp_counts_rep2,supp_counts_rep3):

    """
        This is the extension of the above function that takes the dictinary of interction with their supporing
        pairs and retuns the counts per million of each interactor

    """
    counts_rep1 = supp_count_rep1
    counts_rep2 = supp_counts_rep2
    counts_rep3 = supp_counts_rep3
    final_list ={}

    for k in results_interaction.keys():
        Enh_rep1 = 0
        Enh_rep2 = 0
        Enh_rep3 = 0
        tot_len_Enh = 0

        for values in results_interaction[k]:
            Enh_rep1 = int (values[0]) + Enh_rep1
            Enh_rep2 = int (values[1]) + Enh_rep2
            Enh_rep3 = int (values[2]) + Enh_rep2
            tot_len_Enh = int(values[3])+ tot_len_Enh

        tot_len_Enh = tot_len_Enh
        Enh_rep4 = round(Enh_rep1/(counts_rep1/1000000),3)
        Enh_rep5 = round(Enh_rep2/(counts_rep2/1000000),3)
        Enh_rep6 = round(Enh_rep3/(counts_rep3/1000000),3)
        Enh_rep7 = round((Enh_rep1 * 1000000000) / (counts_rep1 * tot_len_Enh),3)
        Enh_rep8 = round((Enh_rep2 * 1000000000) / (counts_rep2 * tot_len_Enh),3)
        Enh_rep9 = round((Enh_rep3 * 1000000000) / (counts_rep3 * tot_len_Enh),3)
        Enh_rep10 = round(Enh_rep1/(158964165 / 1000000),3)
        Enh_rep11 = round(Enh_rep2/(158701117 / 1000000),3)
        Enh_rep12 = round(Enh_rep2/(156802507 / 1000000),3)

        final_list[k] = [Enh_rep1,Enh_rep2,Enh_rep3,Enh_rep4,Enh_rep5,Enh_rep6,Enh_rep7,Enh_rep8,Enh_rep9,Enh_rep10,Enh_rep11,Enh_rep12,tot_len_Enh]
    return (final_list)

def Main():
    parser = argparse.ArgumentParser(description="General command to manipulate interaction and expresstion data")
    parser.add_argument("PD", help = "PD dataset")
    parser.add_argument("exprs", help ="genes expressed in all celltype")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()

    with open (args.exprs, "r") as expression_file:
        common_genes = expr_gene(expression_file)
    #print(common_genes)

    with open (args.PD, "r") as int_file_PD:
        next(int_file_PD)
        celltype_SP,replicate1_SP_count_PD,replicate2_SP_count_PD  = int_dict_SP(int_file_PD, common_genes)
        celltype_CPM_PD = SP_CPM(celltype_SP, replicate1_SP_count_PD, replicate2_SP_count_PD )

    print(replicate1_SP_count_PD,replicate2_SP_count_PD)


    if args.output:
        with open (args.output, "w") as out:
            for keys,values in celltype_CPM_PD.items():
                out.write (keys+"\t"+str(len(celltype_SP[keys]))+"\t"+"\t".join(str(i) for i in values))
                out.write ("\n")


if __name__ == "__main__":
    Main()
