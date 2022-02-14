# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 18:24:09 2020

@author: Fatemeh
"""

from Bio import SeqIO
import os, sys, re, glob, itertools
import pandas as pd
from operator import itemgetter 
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp



def extract_indecies(list1, list2, list3):
    """extracts the indcies of ids based on their index in another list
    retun the intersection with list3"""
    df_1 = pd.read_csv(list1)
    df_1 = df_1.set_index(df_1.columns[0])
    list_1 = list(df_1.to_numpy().flatten())
    
    df_2 = pd.read_csv(list2)
    df_2 = df_2.set_index(df_2.columns[0])
    list_2 = list(df_2.to_numpy().flatten())
    
    df_3 = pd.read_csv(list3)
    df_3 = df_3.set_index(df_3.columns[0])
    list_3 = list(map(int, list(df_3.to_numpy().flatten())))
    indecies = [i for i, e in enumerate(list_1) if e in set(list_2)]
    lists_intersection = set(indecies).intersection(list_3)
    pd.DataFrame(lists_intersection).to_csv("seq_indecies_to_keep.csv")


def extract_sequences(list_file, fastafilename):
    """Extract the sequences that their index is in the list and save it in a new fasta file"""
    df = pd.read_csv(list_file)
    df = df.set_index(df.columns[0])
    seq_list = list(map(int, list(df.to_numpy().flatten())))
    handle = open(str(list_file)+"_extract_seqs.fasta","w")
    records = list(SeqIO.parse(fastafilename, "fasta"))
    modified_seq = list(itemgetter(*seq_list)(records)) 
    SeqIO.write(modified_seq,handle,"fasta")
    handle.close()
    
    
def extract_k_sequences(k, fastafilename):
    """extract the fist k sequences out of a fasta file"""
    handle = open(str(fastafilename)+"_extract_"+str(k)+"_seqs.fasta","w")
    records = list(SeqIO.parse(fastafilename, "fasta"))
    SeqIO.write(records[:k],handle,"fasta")
    handle.close()
    
    
def remove_seq_with_gap(gap_size, fastafilename):
    """removes the sequences with gaps more than the gap_size"""
    sequence_set = []
    for record in list(SeqIO.parse(fastafilename, "fasta")):
        searchObj = re.finditer(r"[^ATGC]", str(record.seq))  
        i = 0 
        for m in searchObj:
            i = i+1
        if (i < gap_size): #remove seq with more than gap size
            sequence_set.append(record)
    handle = open(str(fastafilename)+"_extract_no_gap_seq.fasta","w")
    SeqIO.write(sequence_set,handle,"fasta")
    handle.close()
    
    
def index_of_seq_with_nogap(gap_size, fastafilename):
    """removes the sequences with gaps more than the gap_size and saves the ids"""
    sequence_set = []
    records = list(SeqIO.parse(fastafilename, "fasta"))
    record_ids = [record.id for record in records]
    for record in records:
        searchObj = re.finditer(r"[^ATGC]", str(record.seq))  
        i = 0 
        for m in searchObj:
            i = i+1
        if (i < gap_size): #remove seq with more than gap size
            sequence_set.append(record_ids.index(record.id))
    pd.DataFrame(sequence_set).to_csv(fastafilename+"_index_no_gap.csv")
            
    
def extract_based_str_epi_index(idx_file, fastafilename): 
    listOfIndecis = pd.read_csv(idx_file)
    listOfIndecis = listOfIndecis['Accession ID']
    sequence_set = []
    handle = open(str(fastafilename).split(".fas")[0] + "_extracted_seqs.fasta","w")
    for record in SeqIO.parse(fastafilename, "fasta"):
        recordStr = record.id
        print(recordStr)
        try:
            EPI_id = (recordStr.split('EPI'))[1].split('|')[0]
            if 'EPI' + EPI_id in listOfIndecis.values:
                sequence_set.append(record)
        except:
            continue
    SeqIO.write(sequence_set, handle,"fasta")
    handle.close()
    return sequence_set


def extract_based_str_virus_name(df, fastafilename, directory): 
    """extract all sequences that their id is in idx_file from fasta fils"""
    # df = pd.read_csv(idx_file)
    df = df['Virus name']
    print(df)
    sequence_set = []
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


def extract_based_country(fastafilename, country): 
    """extract all sequences for a specific country"""
    handle = open(str(fastafilename).split('.fasta')[0] + "_" + country + ".fasta","w")
    for record in SeqIO.parse(fastafilename, "fasta"):
        try:
            str(record.id + record.description).index(country)
            print(record.id + record.description)
            SeqIO.write(record, handle,"fasta")
        except:
            continue
    handle.close()


def extract_cluster_c_seq(labels_file, fasta_files_path, c):
    """extracts the sequences belonging to cluster c (having label c)
    and save it as fasta files"""
    df = pd.read_csv(labels_file)
    df = df.set_index(df.columns[0])
    seq_indecis = df[df['label'] == c].index.values
    for file in glob.glob(os.path.join(fasta_files_path, '*_1004.fasta')):
        extract_based_str_epi_index(seq_indecis, file, c)
        
        
def extract_based_submission_date(start_date, metadata, id_ = 'Virus name'):
    """extract all epi ids of sequences with date after the start_date"""
    metadata_info = pd.read_csv(metadata)
    # tsv_read = pd.read_csv(metadata, sep='\t')
    # metadata_info = tsv_read[[id_, 'Accession ID', 'Submission date']]
    metadata_info = metadata_info.set_index(id_)
    labels = []
    for date in metadata_info['Collection date'].unique():
        try:
            datetime = pd.to_datetime(date)
            if datetime > start_date:
                labels.append(metadata_info.index[metadata_info['Collection date'] == date])
        except:
            continue
    return list(itertools.chain(*labels)), metadata_info
    

def get_strain_index_fasta(fastafilename): 
    index_set = []
    for record in SeqIO.parse(fastafilename, "fasta"):
        index_set.append(record.id)
    pd.DataFrame(index_set).to_csv("/alina-data0/Fatemeh/msa_0511_strain.csv")
    return index_set


def extract_country_ids(csv_file):
    df = pd.read_csv(csv_file)
    countries = df[df.columns[-1]].unique()
    # countries = list(set([x.strip() for x in countries]))
    for country in countries:
        group = df.groupby('1').get_group(country)
        group.to_csv(str(csv_file).split(".csv")[0] + '_' + str(country) + '.csv')
    
    
def replace_space_in_ids(fasta_file):
    """replace the spaces in the sequences id with _"""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    handle = open(str(fasta_file).split('.fasta')[0] + "_space_fixed.fasta","w")
    for record in records:
        my_seq = SeqRecord(record.seq, id = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def groupby_m_y(metadata, tuple_, directory):
    location = tuple_[0].split('_')
    df = tuple_[1]
    dates = df['Collection date']
    df['Collection date'] = pd.to_datetime(dates).dt.to_period('m')
    groupby = df.groupby(by="Collection date")
    counts = groupby.count()['Virus name']
    groups = list(groupby)
    path = directory + location[0] + '/' + location[1] + '/'
    if not os.path.exists(path):
        os.makedirs(path)
    for group in groups:
        file_name = path + tuple_[0] + '_' + group[0].strftime('%Y-%m')
        group[1].to_csv(file_name + '.csv')
        extract_based_str_virus_name(group[1], fasta_file, file_name)
    return counts, groups
    
    
def get_m_y_counts(metadata):
    tsv_read = pd.read_csv(metadata, sep='\t')
    tsv_read = tsv_read[['Accession ID', 'Collection date', 'Virus name', 'Location', 'Host']]
    tsv_read = tsv_read.loc[tsv_read['Host'] == 'Human']
    location_list = []
    for location in tsv_read['Location']:
        continent = location.split('/')[0].strip()
        country = location.split('/')[1].strip()
        # state = location.split('/')[1:].strip()
        location_list.append(continent + '_'+ country)
    tsv_read['Location'] = location_list
    list_ = list(tsv_read.groupby(by="Location"))
    count_list = []
    columns = []
    for group in list_:
        columns.append(group[0])
        # extract_based_str_virus_name(group[1], fasta_file, fasta_file.split('.fasta')[0] + group[0])
        counts, groupby_date = groupby_m_y(metadata, group, '')
        count_list.append(counts)
    count_df = pd.concat(count_list, axis = 1)
    count_df.columns = columns
    count_df.to_csv(metadata.split('.tsv')[0] + 'counts.csv')
    
    
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
    metadata = '/Users/fatemehmohebbi/python_run/continents_metadata/Africa_metadata_1126.csv'
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
        file_name = path + location + '_' + group[0].strftime('%Y-%m')
        # group[1].to_csv(file_name + '.csv')
        extract_based_str_virus_name(group[1], fasta_file, file_name)
    return counts, groups

def parallel(file_path):
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
    print(df)
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
        
        
def seperate_states(fasta_file, metadata):
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
        
def extract_singe_state(fasta_file, metadata,state_):
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
    """continent metadata and path to fasta files for countries of that continent"""
    file_path = sys.argv[1] 
    extract_k_sequences(100, file_path)
    # metadata = sys.argv[2] 
    # df = pd.read_csv(metadata)
    # extract_based_str_virus_name(df,file_path, file_path.split('.fasta')[0] + 'canada')
    # list_ = extract_based_submission_date(pd.to_datetime('12/1/2020'), metadata) #mm/dd/yyyy
    # list_[0].to_csv(str(metadata).split(".csv")[0] + '12_1_2020.csv')
    
