import io 
import os
import sys
from openpyxl.workbook.workbook import Workbook
import pandas as pd
import numpy as np

if not len(sys.argv) ==3:
    print ("\nError:\tincorrect number of command-line arguments")
    print ("Syntax:\tmatrix.py [Input VCF] [Output xlsx]\n")
    sys.exit()

if sys.argv[1] == sys.argv[2]:
    #boolOverwrite = raw_input("\nInput file is the same as the output file - overwrite input file (sim/nao)?\n")
    print ("Error:\tInput file is the same as the output file - choose a different output file\n")
    sys.exit()

input_path = sys.argv[1]
output_path = sys.argv[2]

# input_path = "C:/Users/markc/OneDrive - University of Cambridge/Gupta Lab/HCW Clinical and Job details/sequences/bash testing zone/Hospital_1_vaccinated_rename_split.ann.vcf"
# output_path = "C:/Users/markc/OneDrive - University of Cambridge/Gupta Lab/HCW Clinical and Job details/sequences/bash testing zone/1matrix.xlsx"

def read_vcf(path):
    with open(path,'r') as f:
        lines = [l for l in f if not l.startswith(('\"##','##'))]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    f.close()

ann_vcf = read_vcf(input_path)
#keep first annotation only, as the others are mostly upstream/downstream/intergenic
ann_vcf['INFO'] = ann_vcf['INFO'].str.split(',',expand=True)[0]

first_annotation = ann_vcf['INFO'].str.split('|',expand = True)
raw_annotation_header = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
annotation_header = raw_annotation_header.split(' | ')
first_annotation.columns = annotation_header
first_annotation
#re-insert annotation, gene name, HGVS.c and HGVS.p into the vcf dataframe.
ann_vcf = ann_vcf.merge(first_annotation[['Annotation','Gene_Name','HGVS.c','HGVS.p']],left_index = True, right_index=True)
ann_vcf['AA_change'] = ann_vcf['HGVS.p'].str.replace('p\.','')

#add amino acid change (Xaa) for unknown nucleotide change
ann_vcf['HGVS.p'] = ann_vcf.apply(
    lambda row: ann_vcf.loc[ann_vcf['POS']==row['POS'],'HGVS.p'][ann_vcf.loc[ann_vcf['POS'] == row['POS'],'HGVS.p'].notna()].to_list()[0][:-3] + 'Xaa' if row['HGVS.p'] == None else row['HGVS.p'],
    axis =1  
)
#fill gene_name for unknown nucleotide change
ann_vcf['Gene_Name'] = ann_vcf.apply(
    lambda row: ann_vcf.loc[ann_vcf['POS']==row['POS'],'Gene_Name'][ann_vcf.loc[ann_vcf['POS']==row['POS'],'Gene_Name'].notna()].to_list()[0] if row['Gene_Name'] == None else row['Gene_Name'], 
    axis =1
)

#convert amino acid codes from 3seq to 1seq
#Bio is bugged where PyLance thinks it is non-existent but can still be ran in command line
from Bio import SeqUtils
ann_vcf['From_AA'] = ann_vcf['HGVS.p'].str.replace('p\.','').str.slice(stop=3).apply(SeqUtils.seq1)
ann_vcf['location'] = ann_vcf['HGVS.p'].str.replace('p\.','').str.slice(start=3,stop=-3)
ann_vcf['To_AA'] = ann_vcf['HGVS.p'].str.replace('p\.','').str.slice(start=-3).apply(SeqUtils.seq1)
ann_vcf['AA_mutation'] = ann_vcf['From_AA'] + ann_vcf['location'].astype(str) + ann_vcf['To_AA']
ann_vcf['mutation'] = ann_vcf['REF'] + ann_vcf['POS'].astype(str) + ann_vcf['ALT'] + ',' + ann_vcf['Gene_Name'] +',' + ann_vcf['AA_mutation']

# Fully converted ann_vcf file to matrixable format

v2v_patient_list = ann_vcf.columns[ann_vcf.columns.str.startswith(('Hospital_','hospital_'))]
patient_list = v2v_patient_list.insert(0,"gi|2050056574|gb|MZ359841.1|") # insert delta 
patient_list = patient_list.insert(0,"MN908947.3") # add Wuhan-1 as reference 


def define_direction(row): # check for any * mutations, 2 = REF contains *, 3 = ALT contains *, 4 = mutation to mutation
    if row['diff'] == 1:
        val = '(Forward)'
    elif row['diff'] == -1:
        val = '(Reverse)'
    elif row['diff'] == 2:
        val = '(n_to_mutated)'
    elif row['diff'] ==3:
        val = '(mutated_to_n)'
    elif row['diff']==4:
        val = '(mutated_to_mutated)'
    return val
#fill box
def find_n_included_mutations(patient_a, patient_b):
    mutation_count='placeholder' #to prevent local variable referenced before assignment error
    if patient_a == patient_b:
        output = 'mutation count:0'
        mutation_count = 0
    else:
        diff_a2b = ann_vcf[patient_b] - ann_vcf[patient_a]

        if not any(diff_a2b): # not any() returns true when all diff_a2b values == 0
            output = 'mutation count:0'
            mutation_count = 0
        else: # we made sure we aren't comparing the same patient, or two patients with equal sequences, now resolve differing variants 
            mutations_a2b = ann_vcf.loc[diff_a2b!=0, ['REF','POS','ALT','Gene_Name','From_AA','location','To_AA','mutation','AA_mutation']] # define comparing terms
            mutations_a2b['diff'] = diff_a2b

            dup_variant = mutations_a2b[(mutations_a2b['POS'].duplicated(keep=False))]
            #resolve different variants in the two patients
            for conflicts in range(len(dup_variant['POS'].unique())):
                dup_variant = mutations_a2b[(mutations_a2b['POS'].duplicated(keep=False))]
                resolved_dup_variant = pd.DataFrame.from_dict(dup_variant.to_dict(),orient='index')
                #-1 is from, 1 is To at 'diff'
                from_dup_variant_index = dup_variant.loc[dup_variant['diff']== -1].index[0]
                to_dup_variant_index = dup_variant.loc[dup_variant['diff']==1].index[0]
                resolver_row_index = max(from_dup_variant_index,to_dup_variant_index) + 1

                resolved_dup_variant[resolver_row_index] = resolved_dup_variant.loc[:,to_dup_variant_index]
                resolved_dup_variant.loc['REF',resolver_row_index] = resolved_dup_variant.loc['ALT',from_dup_variant_index]
                resolved_dup_variant.loc['From_AA',resolver_row_index] = resolved_dup_variant.loc['To_AA',from_dup_variant_index]

                resolved_dup_variant.loc['AA_mutation',resolver_row_index] = resolved_dup_variant.loc['From_AA',resolver_row_index] + resolved_dup_variant.loc['location',resolver_row_index] + resolved_dup_variant.loc['To_AA',resolver_row_index]
                resolved_dup_variant.loc['mutation',resolver_row_index] = resolved_dup_variant.loc['REF',resolver_row_index] + str(resolved_dup_variant.loc['POS',resolver_row_index]) + resolved_dup_variant.loc['ALT',resolver_row_index] + ',' + resolved_dup_variant.loc['Gene_Name',resolver_row_index] +',' + resolved_dup_variant.loc['AA_mutation',resolver_row_index]
                resolved_dup_variant.loc['diff',resolver_row_index] = 2 if resolved_dup_variant[resolver_row_index]['REF'] == '*'  else resolved_dup_variant.loc['diff',resolver_row_index]  # check for any * mutations, 2 = REF contains *, 3 = ALT contains *, 4 = mutation to mutation
                resolved_dup_variant.loc['diff',resolver_row_index] = 3 if resolved_dup_variant[resolver_row_index]['ALT'] == '*' else resolved_dup_variant.loc['diff',resolver_row_index]
                resolved_dup_variant.loc['diff',resolver_row_index] = 4 if not resolved_dup_variant.loc[['REF','ALT'],resolver_row_index].str.contains('\*').any() else resolved_dup_variant.loc['diff',resolver_row_index]

                mutations_a2b = mutations_a2b.append(resolved_dup_variant[resolver_row_index]).sort_index().drop([from_dup_variant_index,to_dup_variant_index])

            mutations_a2b['direction'] = mutations_a2b.apply(define_direction,axis=1) # apply function on row 
            mutations_a2b['output_string'] = mutations_a2b['mutation'] + mutations_a2b['direction']
            mutation_count = len(mutations_a2b)
            output = 'mutation_count:' + str(mutation_count) + '\n' + '\n'.join(mutations_a2b['output_string'].apply(str))
    return (output,mutation_count)
#fill row
def n_included_mutation_comparison_against_all(patient_for_comparison): #define function that creates list for filling in rows 
    Mutation_patient_list = []
    mutation_count_patient_list = []
    for patient in patient_list:
        pairwise_mutation = find_n_included_mutations(patient_for_comparison,patient)[0]
        Mutation_patient_list.append(pairwise_mutation)

        pairwise_mutation_count = find_n_included_mutations(patient_for_comparison,patient)[1]
        mutation_count_patient_list.append(pairwise_mutation_count)
    return Mutation_patient_list, mutation_count_patient_list
# ------------------------------------------ DF CREATION ---------------------------------------------------------------#
n_included_df = pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_included_df.loc[anchoring_patient] = n_included_mutation_comparison_against_all(anchoring_patient)[0]

n_included_count_df = pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_included_count_df.loc[anchoring_patient] = n_included_mutation_comparison_against_all(anchoring_patient)[1]
# ------------------------------------------END OF FIRST SET DF CREATION ------------------------------------------------#

# Define functions for exlcuding mutations that involves unknonwn nucleotides

def find_true_mutations(patient_a, patient_b): # excludes unknown nucleotides
    mutation_count='placeholder' #to prevent local variable referenced before assignment error
    if patient_a == patient_b:
        output = 'mutation count:0'
        mutation_count = 0
    else:
        diff_a2b = ann_vcf[patient_b] - ann_vcf[patient_a]

        if not any(diff_a2b): # not any() returns true when all diff_a2b values == 0
            output = 'mutation count:0'
            mutation_count = 0
        else: # we made sure we aren't comparing the same patient, or two patients with equal sequences, now resolve differing variants 
            mutations_a2b = ann_vcf.loc[diff_a2b!=0, ['REF','POS','ALT','Gene_Name','From_AA','location','To_AA','mutation','AA_mutation']] # define comparing terms
            mutations_a2b['diff'] = diff_a2b
            
            dup_variant = mutations_a2b[(mutations_a2b['POS'].duplicated(keep=False))]
            #resolve different variants in the two patients
            for conflicts in range(len(dup_variant['POS'].unique())):
                dup_variant = mutations_a2b[(mutations_a2b['POS'].duplicated(keep=False))]
                resolved_dup_variant = pd.DataFrame.from_dict(dup_variant.to_dict(),orient='index')
                #-1 is from, 1 is To at 'diff'
                from_dup_variant_index = dup_variant.loc[dup_variant['diff']== -1].index[0]
                to_dup_variant_index = dup_variant.loc[dup_variant['diff']==1].index[0]
                resolver_row_index = max(from_dup_variant_index,to_dup_variant_index) + 1

                resolved_dup_variant[resolver_row_index] = resolved_dup_variant.loc[:,to_dup_variant_index]
                resolved_dup_variant.loc['REF',resolver_row_index] = resolved_dup_variant.loc['ALT',from_dup_variant_index]
                resolved_dup_variant.loc['From_AA',resolver_row_index] = resolved_dup_variant.loc['To_AA',from_dup_variant_index]

                resolved_dup_variant.loc['AA_mutation',resolver_row_index] = resolved_dup_variant.loc['From_AA',resolver_row_index] + resolved_dup_variant.loc['location',resolver_row_index] + resolved_dup_variant.loc['To_AA',resolver_row_index]
                resolved_dup_variant.loc['mutation',resolver_row_index] = resolved_dup_variant.loc['REF',resolver_row_index] + str(resolved_dup_variant.loc['POS',resolver_row_index]) + resolved_dup_variant.loc['ALT',resolver_row_index] + ',' + resolved_dup_variant.loc['Gene_Name',resolver_row_index] +',' + resolved_dup_variant.loc['AA_mutation',resolver_row_index]
                resolved_dup_variant.loc['diff',resolver_row_index] = 2 if resolved_dup_variant[resolver_row_index]['REF'] == '*'  else resolved_dup_variant.loc['diff',resolver_row_index]  # check for any * mutations, 2 = REF contains *, 3 = ALT contains *, 4 = mutation to mutation
                resolved_dup_variant.loc['diff',resolver_row_index] = 3 if resolved_dup_variant[resolver_row_index]['ALT'] == '*' else resolved_dup_variant.loc['diff',resolver_row_index]
                resolved_dup_variant.loc['diff',resolver_row_index] = 4 if not resolved_dup_variant.loc[['REF','ALT'],resolver_row_index].str.contains('\*').any() else resolved_dup_variant.loc['diff',resolver_row_index]

                mutations_a2b = mutations_a2b.append(resolved_dup_variant[resolver_row_index]).sort_index().drop([from_dup_variant_index,to_dup_variant_index])
            
            mutations_a2b['direction'] = mutations_a2b.apply(define_direction,axis=1)
            mutations_a2b['output_string'] = mutations_a2b['mutation'] + mutations_a2b['direction']
            mutations_a2b = mutations_a2b[~mutations_a2b.output_string.str.contains("\*")] ### exclude rows where there is a * included
            mutation_count = len(mutations_a2b)
            output = 'mutation_count:' + str(mutation_count) + '\n' + '\n'.join(mutations_a2b['output_string'].apply(str))
    return (output,mutation_count)

def true_mutation_comaprison_against_all(patient_for_comparison): #define function that creates list for filling in rows 
    Mutation_patient_list = []
    mutation_count_patient_list = []
    for patient in patient_list:
        pairwise_mutation = find_true_mutations(patient_for_comparison,patient)[0]
        Mutation_patient_list.append(pairwise_mutation)

        pairwise_mutation_count = find_true_mutations(patient_for_comparison,patient)[1]
        mutation_count_patient_list.append(pairwise_mutation_count)
    return Mutation_patient_list, mutation_count_patient_list

n_not_included_df =  pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_not_included_df.loc[anchoring_patient] = true_mutation_comaprison_against_all(anchoring_patient)[0]

n_not_included_count_df = pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_not_included_count_df.loc[anchoring_patient] = true_mutation_comaprison_against_all(anchoring_patient)[1]

#---------------------------------- write to excel file --------------------------------------------------------------
from openpyxl import Workbook
from openpyxl.reader.excel import load_workbook

# # sheets can not store sheet name with / in it
# from openpyxl.workbook.child import INVALID_TITLE_REGEX
# import re
# title = re.sub(INVALID_TITLE_REGEX,'_',output_path)

with pd.ExcelWriter(output_path, engine = 'openpyxl') as writer:
    #n_included
    n_included_df.to_excel(writer, sheet_name = 'n_included', index = True)
    n_included_count_df.to_excel(writer,sheet_name='n_included_dists',index=True)
    #n_not_included
    n_not_included_df.to_excel(writer, sheet_name='n_not_included', index = True)
    n_not_included_count_df.to_excel(writer,sheet_name='n_not_included_dists',index=True)
    writer.save()

print(f'completed:{output_path}')