# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 21:17:39 2020

@author: Fatemeh
"""

import os, sys, re, glob
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord

# genes 1-based
genes = {
    "S": (24035, 28971 + 1, "S"),
    #"S": (23949+1, 29303+1, "S"), for msa_1126 
    # "S": (21563, 25384, "S"),
    "ORF_3a": (25393, 26220, "ORF_3a"),
    "E": (26245, 26472, "E"),
    "M": (26523, 27191, "M"),
    "ORF_6": (27202, 27387, "ORF_6"),
    "ORF": (27394, 27759, "ORF"),
    "ORF_8": (27894, 28259, "ORF_8"),
    "N": (28274, 29533, "N")
}


def get_gene(seq, gene):
    return seq[genes[gene][0] - 1: genes[gene][1] - 1]


def codon_to_amino(codon):
    if "-" in codon or len(codon) < 3:
        return '_'
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    return table[codon]


def extract_EPI_id(labels):
    ids = []
    for label in labels:
        try:
            ids.append(label[label.find("EPI"):].split('|')[0])
        except:
            ids.append(label)
    ids = np.asarray(ids)
    return ids


def extract_strain_name(labels):
    ids = []
    for label in labels:
        try:
            ids.append(label.split("hCoV-19/")[-1].split('|')[0])
        except:
            ids.append(label)
    ids = np.asarray(ids)
    return ids


def extract_collection_date(labels):
    dates = []
    for label in labels:
        try:
            date_str = label.split('|')[-2]
            if len(re.findall(r'[^0-9]', date_str)) == 2:
                valid = pd.to_datetime(date_str)
                dates.append(date_str)
        except:
            dates.append(np.nan)
    dates = np.asarray(dates)
    return dates


def translate(seq):
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:
            break
        protein += codon_to_amino(codon)
    return protein


def find_mutation_alignment(fastaFile, refrenceFile):
    """aligns all the sequences in fastafile to the refrence and looks for the mutation D614G"""
    ref = list(SeqIO.parse(refrenceFile, "fasta"))
    records = list(SeqIO.parse(fastaFile, "fasta"))
    mutation_state_arr = np.zeros(len(records))
    counter = 0
    ids = []
    for record in records:
        print(counter)
        ids.append(record.id)
        alignments = pairwise2.align.localms(ref[0].seq, record.seq, 1, -2, -2, -2, one_alignment_only = 'True')
        alignment_str = (format_alignment(*alignments[0],full_sequences=True))
        [target, match, query, score, blank] = str(alignment_str).split('\n')
        hapl_s_gene = get_gene(Seq(query),'S')
        try:
            hapl_s_amino = translate(hapl_s_gene)
            if (hapl_s_amino[613] == 'G'):
                mutation_state_arr[counter] = 1
                print("mutation happened")
        except:
            continue
        counter = counter + 1
    print(ids)
    print(mutation_state_arr)
    pd.DataFrame({'D614G_mutation': mutation_state_arr, 'id': ids}).to_csv(str(fastaFile).split(".fas")[0]+"_mutation_st.csv")
 
    
def find_D614G_mutation(fastaFile):
    """looks for the mutation D614G in all the sequences in the fastafile. First it extracts the S gene
    Then translate it to protein and look for D614G"""
    records = list(SeqIO.parse(fastaFile, "fasta"))
    mutation_state_arr = np.zeros(len(records))
    counter = 0
    ids = []
    for record in records:
        print(counter)
        ids.append(record.id + record.description)
        hapl_s_gene = get_gene(record.seq,'S')
        try:
            hapl_s_amino = translate(hapl_s_gene)
            if (hapl_s_amino[613] == 'G'):
                mutation_state_arr[counter] = 1
                print("mutation happened in: ", record.id)
                print(hapl_s_amino)
        except:
            continue
        counter = counter + 1
    ids = extract_strain_name(ids)
    pd.DataFrame({'D614G_mutation': mutation_state_arr, 'id': ids}).to_csv(str(fastaFile).split(".fas")[0]+"_mutation_st.csv")
 
    
    
def find_Delta(fastaFile):
    """looks for the mutation D614G in all the sequences in the fastafile"""
    records = list(SeqIO.parse(fastaFile, "fasta"))
    mutation_state_arr = np.zeros(len(records))
    handle = open(str(fastaFile).split(".fasta")[0]+"_delta_seq.fasta","w")
    counter = 0
    ids = []
    for record in records:
        print(counter)
        ids.append(record.id + record.description)
        try:
            hapl_s_amino = record.seq
            if (hapl_s_amino[613] == 'G') & (hapl_s_amino[451] == 'R') & \
                (hapl_s_amino[477] == 'K') & (hapl_s_amino[680] == 'R') & \
                 (hapl_s_amino[949] == 'N') & (hapl_s_amino[18] == 'R') & \
                     (hapl_s_amino[141] == 'D') & (hapl_s_amino[94] == 'I'):
                   # (hapl_s_amino[155] == '_') & (hapl_s_amino[156] == '_') & (hapl_s_amino[157] == 'G') & \
                mutation_state_arr[counter] = 1
                print("mutation happened in: ", record.id)
                SeqIO.write(record, handle, "fasta") 
        except:
            continue
        counter = counter + 1
    handle.close()
    epi_id = extract_EPI_id(ids)
    dates = extract_collection_date(ids)
    ids = extract_strain_name(ids)
    pd.DataFrame({'mutation': mutation_state_arr, 'Collection date': dates, \
                  'Accession id': epi_id, 'Virus name': ids}).to_csv(str(fastaFile).split(".fas")[0]+"_B16172.csv") 
    
 
    
def find_pos_in_gene(fasta_file, gene):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        hapl_s_gene = get_gene(record.seq, gene)
        try:
            hapl_s_amino = translate(hapl_s_gene)
            print(hapl_s_amino)
            print('________________________________________')
        except:
            continue
 

def extract_s_gene(fastaFile):
    """extracts s gene from all the sequences in fasta file"""
    counter = 0
    handle = open(str(fastaFile).split(".fasta")[0]+"_s_gene_protein.fasta","w")
    for record in SeqIO.parse(fastaFile, "fasta"):
        print(counter)
        counter = counter + 1
        hapl_s_gene = get_gene(str(record.seq),'S')
        my_seq = SeqRecord(Seq(hapl_s_gene), id = record.id, description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def translate_to_protein(fastaFile):
    """extracts s gene from all the sequences in fasta file"""
    records = list(SeqIO.parse(fastaFile, "fasta"))
    counter = 0
    handle = open(str(fastaFile).split(".fasta")[0]+"_protein.fasta","w")
    for record in records:
        print(counter)
        hapl_s_amino = translate(record.seq)
        counter = counter + 1
        my_seq = SeqRecord(Seq(hapl_s_amino), id = record.id,  description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()


def find_pairs(fastaFile):
    """looks for the mutation D614G in all the sequences in the fastafile"""
    records = list(SeqIO.parse(fastaFile, "fasta"))
    mutation_state_arr = np.zeros(len(records))
    counter = 0
    ids = []
    for record in records:
        print(counter)
        ids.append(record.id + record.description)
        try:
            hapl_s_amino = record.seq
            if (hapl_s_amino[497] != 'R') & (hapl_s_amino[500] != 'N'):
                mutation_state_arr[counter] = 1
                print("mutation happened in: ", record.id)
        except:
            continue
        counter = counter + 1
    epi_id = extract_EPI_id(ids)
    dates = extract_collection_date(ids)
    ids = extract_strain_name(ids)
    pd.DataFrame({'mutation': mutation_state_arr, 'Collection date': dates, \
                  'Accession id': epi_id, 'Virus name': ids}).to_csv(str(fastaFile).split(".fas")[0]+"_498_501_st.csv") 
    
def return_second_part_of_seq(fastaFile, x):
    records = list(SeqIO.parse(fastaFile, "fasta"))
    handle = open(str(fastaFile).split(".fasta")[0]+ str(x) +".fasta","w")
    for record in records:
        my_seq = record[:x]
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
        
if __name__ == '__main__':
    # directory = '/Users/fatemehmohebbi/Downloads/'
    file_path =sys.argv[1]  #directory + 'msa_0207_extract_100_seqs.fasta' 
    # for file in glob.glob(os.path.join(file_path, '*.fasta')):
    extract_s_gene(file_path)
    # fasta_file = sys.argv[1] 
    # fasta_file = "/Users/fatemehmohebbi/Downloads/metadataEuropeBelgium_extract_seqs.fasta"
    # find_mutation(fasta_file, 'N', 203, 'N')