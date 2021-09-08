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


def define_direction(row):
    if row['diff'] == 1:
        val = '(Forward)'
    elif row['diff'] == -1:
        val = '(Reverse)'
    else:
        val = '(what)'
    return val


def find_n_included_mutations(patient_a, patient_b):
    if patient_a == patient_b:
        output = 'mutation count:0'
    else:
        diff_a2b = ann_vcf[patient_b] - ann_vcf[patient_a]
        if not any(diff_a2b):
            output = 'mutation count:0'
        else:
            mutations_a2b = ann_vcf.loc[diff_a2b!=0, ['Gene_Name','mutation']] # define comparing terms
            mutations_a2b['diff'] = diff_a2b
            
            mutations_a2b['direction'] = mutations_a2b.apply(define_direction,axis=1)
            mutations_a2b['output_string'] = mutations_a2b['mutation'] + mutations_a2b['direction']

            mutation_count = len(mutations_a2b)
            output = 'mutation_count:' + str(mutation_count) + '\n' + '\n'.join(mutations_a2b['output_string'].apply(str))
    return output

def n_included_mutation_comparison_against_all(patient_for_comparison): #define function that creates list for filling in rows 
    Mutation_patient_list = []
    for patient in patient_list:
        pairwise_mutation = find_n_included_mutations(patient_for_comparison,patient)
        Mutation_patient_list.append(pairwise_mutation)
    return Mutation_patient_list

n_included_df = pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_included_df.loc[anchoring_patient] = n_included_mutation_comparison_against_all(anchoring_patient)

def find_true_mutations(patient_a, patient_b): # excludes unknown nucleotides
    if patient_a == patient_b:
        output = 'mutation count:0'
    else:
        diff_a2b = ann_vcf[patient_b] - ann_vcf[patient_a]
        if not any(diff_a2b):
            output = 'mutation count:0'
        else:
            mutations_a2b = ann_vcf.loc[diff_a2b!=0, ['Gene_Name','mutation']] # define comparing terms
            mutations_a2b['diff'] = diff_a2b
            
            mutations_a2b['direction'] = mutations_a2b.apply(define_direction,axis=1)
            mutations_a2b['output_string'] = mutations_a2b['mutation'] + mutations_a2b['direction']
            mutations_a2b = mutations_a2b[~mutations_a2b.output_string.str.contains("\*")] # exclude rows where there is a * included

            mutation_count = len(mutations_a2b)
            output = 'mutation_count:' + str(mutation_count) + '\n' + '\n'.join(mutations_a2b['output_string'].apply(str))
    return output

def true_mutation_comaprison_against_all(patient_for_comparison): #define function that creates list for filling in rows 
    Mutation_patient_list = []
    for patient in patient_list:
        pairwise_mutation = find_true_mutations(patient_for_comparison,patient)
        Mutation_patient_list.append(pairwise_mutation)
    return Mutation_patient_list

n_not_included_df =  pd.DataFrame(columns=patient_list,index = patient_list)
for anchoring_patient in patient_list:
    n_not_included_df.loc[anchoring_patient] = true_mutation_comaprison_against_all(anchoring_patient)

#write to excel file
from openpyxl import Workbook
from openpyxl.reader.excel import load_workbook

# # sheets can not store sheet name with / in it
# from openpyxl.workbook.child import INVALID_TITLE_REGEX
# import re
# title = re.sub(INVALID_TITLE_REGEX,'_',output_path)

with pd.ExcelWriter(output_path, engine = 'openpyxl') as writer:
    n_included_df.to_excel(writer, sheet_name = 'n_included', index = True)
    n_not_included_df.to_excel(writer, sheet_name='n_not_included', index = True)
    writer.save()