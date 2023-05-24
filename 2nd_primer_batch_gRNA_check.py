# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:25:51 2022

Design primers for group of sequences such that the gRNA target sequence is in the amplicon.

Input:
    1. CSV file with gRNA names and their target sequences 
    NOTE: The extra sequences in front of the sequences need to be removed - the program should be altered accordingly
    2. Parameters for primer design should be specified.
    3. Input folder with the target genes for each gRNA in a separate fasta file named accrdingly - e.g: A1.fasta with 1 or more sequence
    
Output:
    1. Test if gRNA target is in sequence
    2. Design primers to amplify the target sequences in the input genes

@author: Amrita
"""
from Bio import SeqIO
from os import listdir
import primer3
import csv

#Read all input data

input_folder = "D:/Projects/NUS - CRISPR/Primer design/2nd PCR/V1/genes/"
input_csv = "D:/Projects/NUS - CRISPR/Primer design/2nd PCR/V1/NLR_ploop_KO_oligoorder.csv"

input_filenames = listdir(input_folder)

csv_data = []
with open(input_csv) as csvfile:
    filereader = csv.reader(csvfile)
    for row in filereader:
        if row[0]+'.fasta' in input_filenames:
            csv_data.append(row)
            
        #Exclude plate 2,3 from analysis!
        if 'Plate 2' in row:
            break
        
csvfile.close()


#Design primers for all genes
output_primers = {}
targeted_sequences = {}

flanking_region = 100 #Change according to extra region around target sequence that should be amplified

for f in input_filenames:    
    input_sequences = list(SeqIO.parse(open(input_folder+f), "fasta"))
    grna = f.split('.')[0]
    
    output_primers[grna] = []
    
    for template_sequence in input_sequences:
        
        #Identify region containing target sequence
                
        for row in csv_data:
            if row[0] == grna:
                target_sequence = row[2][5:]
                start_target_sequence = str(template_sequence.seq).index(target_sequence)
                
        target_range = [start_target_sequence - flanking_region, len(target_sequence) + 2*flanking_region]
        
        
        return_primers = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': str(template_sequence.id),
            'SEQUENCE_TEMPLATE': str(template_sequence.seq),
            'SEQUENCE_INCLUDED_REGION': target_range
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 15,
            'PRIMER_MAX_SIZE': 26,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[200, 270]],
        })
        
        output_primers[grna].append([return_primers, template_sequence.id])
        
        targeted_sequences[str(template_sequence.id)] = [str(template_sequence.seq)[(start_target_sequence - 5):start_target_sequence+1], str(template_sequence.seq)[len(target_sequence) + 2*flanking_region + 1: len(target_sequence) + 2*flanking_region + 6]]
        
final_primer_set = []

for key in output_primers:
    
    for primers in output_primers[key]:
    
        left_primer = primers[0]['PRIMER_LEFT_0_SEQUENCE']
        right_primer = primers[0]['PRIMER_RIGHT_0_SEQUENCE']
        left_temp = primers[0]['PRIMER_LEFT_0_TM']
        right_temp = primers[0]['PRIMER_RIGHT_0_TM']
    
        final_primer_set.append([key+'_'+primers[1], left_primer, left_temp, right_primer, right_temp])
    
    
#Write to output csv file

# output_file = "D:/Projects/NUS - CRISPR/Primer design/2nd PCR/V1/2nd_pcr_primers.csv"
# csv_columns = ['gRNA', 'Left primer', 'Tm', 'Right primer', 'Tm']

# with open(output_file, 'w') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(csv_columns)
    
#     for row in final_primer_set:
#         writer.writerow(row)
        
# csvfile.close()
    
    
    
    
