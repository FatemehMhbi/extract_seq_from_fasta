#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 11:39:47 2021

@author: fatemehmohebbi
"""

from Bio import SeqIO
import os, sys, re, glob, itertools
import pandas as pd
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import random
from extract_sequences_by_id_k_fasta import extract_based_str_virus_name, extract_based_country, \
    get_epi_id_fasta


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

    
def groupby_m_y(tuple_, directory):
    location = tuple_[0].split('_')
    df = tuple_[1]
    fasta_ids = get_epi_id_fasta(fasta_file)
    ixs = set(df['Accession ID']).intersection(fasta_ids['Accession ID'])
    df.index = df['Accession ID']
    df = df.loc[ixs]
    dates = df['Collection date']
    df['Collection date'] = pd.to_datetime(dates).dt.to_period('m')
    groupby = df.groupby(by="Collection date")
    counts = groupby.count()[df.columns[0]]
    groups = list(groupby)
    path = directory + location[0] + '/'
    if not os.path.exists(path):
        os.makedirs(path)
    for group in groups:
        file_name = path + tuple_[0] + '_' + group[0].strftime('%Y-%m')
        # group[1].to_csv(file_name + '.csv')
        random_samples = random_pick(group[1], 1000, fasta_file)
        handle =  open(file_name + ".fasta","w")
        records = extract_based_str_virus_name(random_samples, fasta_file, file_name)
        SeqIO.write(records, handle,"fasta")
        handle.close()
    return counts, groups
    
    
def get_m_y_counts(metadata, region):
    """region is an int, 0 for continent, 1 for country, 2 for state"""
    tsv_read = pd.read_csv(metadata) #, sep='\t')
    tsv_read = tsv_read[['Accession ID', 'Collection date', 'Virus name', 'Location', 'Host']]
    tsv_read = tsv_read.drop_duplicates(subset=['Accession ID'])
    tsv_read = tsv_read.loc[tsv_read['Host'] == 'Human']
    tsv_read = remove_no_val_dates(tsv_read)
    location_list = []
    for location in tsv_read['Location']:
        region_ = location.split('/')[region].strip()
        location_list.append(region_) #(continent + '_'+ country)
    tsv_read['Location'] = location_list
    list_ = list(tsv_read.groupby(by="Location"))
    count_list = []
    columns = []
    for group in list_:
        columns.append(group[0])
        # extract_based_str_virus_name(group[1], fasta_file, fasta_file.split('.fasta')[0] + group[0])
        counts, groupby_date = groupby_m_y(group, '')
        count_list.append(counts)
    count_df = pd.concat(count_list, axis = 1)
    count_df.columns = columns
    count_df.to_csv(metadata.split('.csv')[0] + '_m_y_counts.csv')
    
    
def save_fasta_from_id_list(file_path, fasta_file):
    for file in glob.glob(os.path.join(file_path, '*.csv')):
      extract_based_str_virus_name(file, fasta_file)


def seperate_continent_metadata(metadata):
    meta_df = pd.read_csv(metadata)
    location_list = []
    for location in meta_df["Location"]:
        continent = location.split('/')[0].strip()
        location_list.append(continent) # + '_'+ country)
    meta_df['continent'] = location_list
    for continent in meta_df['continent'].unique():
        meta_df.loc[meta_df['continent'] == continent].to_csv(continent + '_metadata_1126.csv')


def seperate_countries_fasta(fasta_file, metadata):
    meta_df = pd.read_csv(metadata)
    location_list = []
    for location in meta_df["Location"]:
        try:
            country = location.split('/')[1].strip()
            location_list.append(country)
        except:
            location_list.append('unkonwn')
            continue
    meta_df['country'] = location_list
    for country in meta_df['country'].unique():
        extract_based_country(fasta_file, country)


def seperat_fasta_m_y(fasta_file):
    location = fasta_file.split('.fasta')[0].split('_')[-1]
    print(location)
    # metadata = '/Users/fatemehmohebbi/python_run/continents_metadata/Africa_metadata_1126.csv'
    df = pd.read_csv(metadata)
    dates = df['Collection date']
    df['Collection date'] = pd.to_datetime(dates).dt.to_period('m')
    groupby = df.groupby(by="Collection date")
    counts = groupby.count()['Virus name']
    groups = list(groupby)
    path = fasta_file[:fasta_file.rfind('/') + 1] + location + '/'
    if not os.path.exists(path):
        os.makedirs(path)
    for group in groups:
        file_name = path + location + '_all_' + group[0].strftime('%Y-%m')
        # group[1].to_csv(file_name + '.csv')
        extract_based_str_virus_name(group[1], fasta_file, file_name)
    return counts, groups


def folder_read_seperate_m_y(file_path):
    """read fasta files in a folder in parallel and save them in seperate fasta file based
    month and year"""
    a_pool = mp.Pool(processes = 5)
    counts = []
    groups = []
    for count, group in a_pool.map(seperat_fasta_m_y, glob.glob(os.path.join(file_path, '*.fasta'))):
        counts.append(count)
        groups.append(group)
        
        
def extract_based_str_virus_name_parallel(group): 
    """extract all sequences that their id is in idx_file from fasta fils"""
    df = group[1]
    fastafilename = '/alina-data0/Fatemeh/Epistasis/MSA/msa_1126/countries/North_America/USA_states/msa_1126_s_gene_ambiguous_replaced_w_blank_America_USA.fasta'
    df = df['Virus name']
    sequence_set = []
    directory = fastafilename.split('.fasta')[0] + group[0] + 'parallel'
    handle =  open(directory + ".fasta","w")
    for record in SeqIO.parse(fastafilename, "fasta"):
        id_ = 'hCoV-19/' + str(record.id + record.description).split('|')[0].split('hCoV-19/')[-1]
        # print(id_)
        if id_ in df.values:
            # sequence_set.append(record)
            SeqIO.write(record, handle,"fasta")
    handle.close()
    if os.path.getsize(directory + ".fasta") == 0:
        os.remove(directory + ".fasta")
    return sequence_set
        
        
def seperate_usa_states(fasta_file, metadata):
    tsv_read = pd.read_csv(metadata)
    tsv_read = tsv_read[['Accession ID', 'Collection date', 'Virus name', 'Location', 'Host']]
    tsv_read = tsv_read.loc[tsv_read['Host'] == 'Human']
    location_list = []
    for location in tsv_read['Location']:
        try:
            state = location.split('/')[2].strip()
            location_list.append(state)
        except:
            location_list.append('unkonwn')
            continue
    tsv_read['state'] = location_list
    list_ = list(tsv_read.groupby(by="state"))
    a_pool = mp.Pool(processes = 10)
    sequences = []
    for seq in a_pool.map(extract_based_str_virus_name_parallel, list_):
        sequences.append(seq)
        
def extract_single_state(fasta_file, metadata,state_):
    tsv_read = pd.read_csv(metadata)
    tsv_read = tsv_read[['Accession ID', 'Collection date', 'Virus name', 'Location', 'Host']]
    tsv_read = tsv_read.loc[tsv_read['Host'] == 'Human']
    location_list = []
    for location in tsv_read['Location']:
        try:
            state = location.split('/')[2].strip()
            location_list.append(state)
        except:
            location_list.append('unkonwn')
            continue
    tsv_read['state'] = location_list
    tsv_read = tsv_read.loc[tsv_read['state'] == state_]
    extract_based_str_virus_name(tsv_read, fasta_file, fasta_file.split('.fasta')[0] + state_)
    

if __name__ == '__main__':
    """metadata and fasta_file should be for same continent or country"""
    fasta_file = sys.argv[1] 
    metadata = sys.argv[2] 
    seperat_fasta_m_y(fasta_file)
    
