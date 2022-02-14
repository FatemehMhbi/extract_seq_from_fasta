#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 11:55:01 2021

@author: fatemehmohebbi
"""

from Bio import SeqIO
import os, sys, re, glob, itertools
import pandas as pd
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import random
import numpy as np
from extract_sequences_by_id_k_fasta import extract_based_str_virus_name, extract_based_country, \
    get_epi_id_fasta


def return_gaps_measure(sequences_file):
    sequence_id = []
    df = []
    for record in sequences_file: #list(SeqIO.parse(fastafilename, "fasta")):
        # id_ = 'hCoV-19/' + str(record.id + record.description).split('|')[0].split('hCoV-19/')[-1]
        length = len(record.seq)
        id_ = 'EPI' + str(record.id + record.description).split('EPI')[1].split('|')[0]
        searchObj = re.finditer(r"[^ATGCactg]", str(record.seq))  
        i = 0 
        for m in searchObj:
            i = i+1
        df.append(length / (i + 1))
        sequence_id.append(id_)
    norm = [float(i)/sum(df) for i in df]
    return pd.DataFrame(norm, index = sequence_id, columns = ['prob'])


def random_pick_probablity(df, n, fasta_file):
    """n is number of the random samples"""
    if df.shape[0] < n:
        samples = df
        print("less than 1000")
    else:
        seq = list(SeqIO.parse(fasta_file, "fasta"))
        probablity = return_gaps_measure(seq)
        ixs = set(df.index).intersection(probablity.index)
        df = df.loc[ixs].dropna()
        df = pd.concat([df, probablity], axis=1)
        random_index = np.random.choice(range(df.shape[0]), n, replace=False, p=df['prob'])
        random.sample(range(df.shape[0]), n)
        samples = df.iloc[random_index]
    return samples


def remove_seq_with_gap(gap_size, sequences_file):
    """removes the sequences with gaps more than the gap_size"""
    sequence_id = []
    df = []
    for record in sequences_file: #list(SeqIO.parse(fastafilename, "fasta")):
        # id_ = 'hCoV-19/' + str(record.id + record.description).split('|')[0].split('hCoV-19/')[-1]
        id_ = 'EPI' + str(record.id + record.description).split('EPI')[1].split('|')[0]
        searchObj = re.finditer(r"[^ATGCactg]", str(record.seq))  
        i = 0 
        for m in searchObj:
            i = i+1
        if (i < gap_size): #remove seq with more than gap size
            df.append(i)
            sequence_id.append(id_)
    return pd.DataFrame(df, index = sequence_id)


def remove_no_val_dates(metadata):
    dates = metadata['Collection date']
    valid_array = []
    for date in dates:
        if len(date) < 6:
            valid_array.append(False)
        else:
            valid_array.append(True)
    return metadata.loc[valid_array]


def random_pick(df, n, fasta_file):
    """n is number of the random samples"""
    if df.shape[0] < n:
        samples = df
        print("less than 1000")
    else:
        seq = extract_based_str_virus_name(df, fasta_file, '')
        df_no_gap = remove_seq_with_gap(60, seq)
        if len(df_no_gap) < n:
            random_index = random.sample(range(df.shape[0]), n)
            samples = df.iloc[random_index]
            print("if removed less than 1000 seq")
        else:
            print("removed")
            random_index = random.sample(range(df.loc[df_no_gap.index].shape[0]), n)
            samples = df.loc[df_no_gap.index].iloc[random_index]
    return samples

    
def groupby_m_y(df, directory, fasta_file):
    random_samples = random_pick_probablity(df, 1000, fasta_file)
    file_name = directory + fasta_file.split('.fasta')[0] + 'samples'
    print(file_name)
    handle =  open(file_name + ".fasta","w")
    records = extract_based_str_virus_name(random_samples, fasta_file, file_name)
    SeqIO.write(records, handle,"fasta")
    handle.close()
    
    
def get_m_y_counts(metadata, fasta_file):
    tsv_read = pd.read_csv(metadata) #, sep='\t')
    tsv_read = tsv_read.drop_duplicates(subset=['Accession ID'])
    tsv_read.index = tsv_read['Accession ID']
    tsv_read = remove_no_val_dates(tsv_read)
    fasta_ids = get_epi_id_fasta(fasta_file)
    ixs = set(tsv_read['Accession ID']).intersection(fasta_ids['Accession ID'])
    tsv_read = tsv_read.loc[ixs]
    groupby_m_y(tsv_read, '', fasta_file)
    # counts.to_csv(metadata.split('.csv')[0] + '_m_y_counts.csv')

    

if __name__ == '__main__':
    """metadata and fasta_file should be for same continent or country"""
    path = sys.argv[1] 
    # metadata = 'metadata_m_y/North America_2021-11.csv' #sys.argv[2] 
    for file in glob.glob(os.path.join(path, '*.fasta')):
        print(file)
        metadata_path = file.split('.fasta')[0].split('_')[-1]
        for metadata in glob.glob(os.path.join(path, '*' + metadata_path + '.csv')):
            print(metadata)
            get_m_y_counts(metadata, file)
            
    
