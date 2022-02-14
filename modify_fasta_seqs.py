# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:34:44 2020

@author: Fatemeh
"""


from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Seq import Seq
import sys, re, os
import numpy as np


def replace_non_actg(string):
    full_pattern = re.compile('[^ACTGactg]')
    return re.sub(full_pattern, '-', string)


def length_of():
    fastafilename = "/Users/fatemehmohebbi/Downloads/ref_s_gene.fasta"
    records = SeqIO.parse(fastafilename, "fasta")
    for record in records:
        print(len(record.seq))
        
        
def take_part_of_seq(fastafilename):
    """extracts between pos and pos2, to extract s gene"""
    records = list(SeqIO.parse(fastafilename, "fasta"))
    length = len(records[0].seq)
    pos = 22000 + 1500 + 500 + 29 + 6
    pos2 = 27000 + 3000 - 500 - 600 + 100 - 29
    print(pos)
    print(pos2)
    handle = open(str(fastafilename).split('.fasta')[0] + "_part.fasta","w")
    my_seq = SeqRecord(Seq(str(records[0].seq[pos : pos2])), id = records[0].id, description = records[0].description)
    SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def join_strings(string, index_list):
    """index_list  is a list of tuples"""
    new_string = ''
    for tuple_ in index_list:
        new_string = new_string + string[tuple_[0]:tuple_[1]]
    return new_string
    
    
def close_gaps_based_ref(fastafilename):
    """remove the gaps based on the gaps on the first sequence"""
    records = list(SeqIO.parse(fastafilename, "fasta"))
    first_seq = records[0].seq
    string_seq = str(first_seq)
    length = len(first_seq)
    print(length)
    string_list = list(filter(None, string_seq.split('-')))
    index_list = []
    start = 0
    end = 0
    for string in string_list:
        start = end + string_seq[end:].find(string)
        end = start + len(string)
        index_list.append((start, end))
    handle = open(str(fastafilename).split('.fasta')[0] + "_gaps_removed.fasta","w")
    for record in records:
        new_seq = join_strings(record.seq, index_list)
        my_seq = SeqRecord(Seq(str(new_seq)), id = record.id, description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
        
    
def close_gaps_based_threshold(fasta_file, threshold):
    print("reading the fasta file ...")
    alignment = AlignIO.read(fasta_file, "fasta")
    num_of_seq = len(alignment)
    print("len(alignment) ", len(alignment))
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    my_pssm = summary_align.pos_specific_score_matrix(consensus)
    print(len(alignment[0]))
    sequence_set = []
    handle = open(str(fasta_file).split(".fasta")[0]+"_gaps_rm_90.fasta","w")
    pos_to_remove = []
    for pos in range(len(alignment[0])):
        print(my_pssm[pos])
        frequency = sum(list(my_pssm[pos].values()))
        if frequency < (num_of_seq * threshold) / 100:
            pos_to_remove.append(pos)
    pos_to_keep = list(set(range(len(alignment[0])))-set(pos_to_remove))
    for align in alignment:
        modified_seq = align[:0]
        for pos in pos_to_keep:
            modified_seq = modified_seq + align[pos]
        sequence_set.append(modified_seq)
    SeqIO.write(sequence_set,handle,"fasta") 
    handle.close()
    
    
def find_intervals(index_list):
    tuples = []
    start = index_list[0] if index_list[0] != 0 else 0
    i = 0
    while i < len(index_list) - 1:
        # print(i)
        if (index_list[i] !=  index_list[i + 1] - 1):
            end = index_list[i] + 1
            tuples.append((start, end))
            start = index_list[i + 1] 
        i = i + 1
        if i == len(index_list) - 1:
            end = index_list[i] + 1
            tuples.append((start, end))  
    return tuples
     
    
def close_gaps_based_threshold_2(seq_file, threshold):
    """removes positions with too many blank"""
    sequences = list(SeqIO.parse(seq_file, "fasta"))
    num_of_seq = len(sequences)
    handle = open(str(seq_file)+"_gap_removed_" + str(threshold) + ".fasta","w")
    summation = np.zeros(len(sequences[0].seq))
    all_blank = "-" * len(sequences[0].seq)
    for seq in sequences:
        matches = np.array([int(i==j) for i,j in zip(seq.seq, all_blank)])
        summation = summation + matches
    print(summation)
    pos_to_remove = np.where(summation < (num_of_seq * threshold) / 100)
    index_list = find_intervals(pos_to_remove[0])
    for record in sequences:
        print(record.id)
        new_seq = join_strings(record.seq[1:], index_list)
        my_seq = SeqRecord(Seq(str(new_seq)), id = record.id, description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def replace_amb_chars_with_blank(seq_file):
    handle = open(str(seq_file).split(".fasta")[0]+"_ambiguous_replaced_w_blank.fasta","w")
    for record in SeqIO.parse(seq_file, "fasta"):
        new_seq = replace_non_actg(str(record.seq))
        my_seq = SeqRecord(Seq(new_seq), id = record.id, description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def change_names(seq_file):
    handle = open(str(seq_file).split(".fasta")[0]+"_ids_changed_.fasta","w")
    count = 1
    for record in SeqIO.parse(seq_file, "fasta"):
        seq_number = str(count) + '_' #record.id.split('_')[-1] + '_'
        # count = count + 1
        # patient_id = record.id.split('_')[0]
        # new_id = patient_id + '_' + seq_number + record.id.split(patient_id)[1]
        new_id = record.id[:10]
        my_seq = SeqRecord(record.seq, id = new_id, description = '')
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def remove_seq_based_name(seq_file):
    """removes all sequences with less than 5 threshold copy number,stated in id"""
    handle = open(str(seq_file).split(".fasta")[0]+"_threshold_2.fasta","w")
    count = 1
    for record in SeqIO.parse(seq_file, "fasta"):
        copy_number = int(record.id.split('_')[-1])
        count = count + 1
        if copy_number >= 2:
            SeqIO.write(record, handle, "fasta") 
    handle.close()   
    
    
def remove_blanks_at_begining(fastafilename):
    """remove blanks at the begining of aligned sequences"""
    records = list(SeqIO.parse(fastafilename, "fasta"))
    length = len(records[0].seq)
    handle = open(str(fastafilename).split('.fasta')[0] + "_half.fasta","w")
    for record in records:
        my_seq = SeqRecord(Seq(str(record.seq[2:])), id = record.id, description = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def print_seq(fastafilename):
    """remove blanks at the begining of aligned sequences"""
    records = list(SeqIO.parse(fastafilename, "fasta"))
    for record in records:
        if record.seq[0] == '-':
            print(record.seq[:10])
            print(record.id)


    
if __name__ == '__main__':
    #directory = '/Users/fatemehmohebbi/Downloads/'
    fastafilename = sys.argv[1] 
    close_gaps_based_threshold_2(fastafilename, 90)
    
    
