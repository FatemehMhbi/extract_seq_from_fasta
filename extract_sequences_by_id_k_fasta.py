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
    df = df['Accession ID']
    sequence_set = []
    # handle =  open(directory + ".fasta","w")
    for record in SeqIO.parse(fastafilename, "fasta"):
        # id_ = 'hCoV-19/' + str(record.id + record.description).split('|')[0].split('hCoV-19/')[-1]
        id_ = 'EPI' + str(record.id + record.description).split('EPI')[1].split('|')[0]
        # print(id_)
        if id_ in df.values:
            sequence_set.append(record)
            # SeqIO.write(record, handle,"fasta")
    # handle.close()
    # if os.path.getsize(directory + ".fasta") == 0:
    #     os.remove(directory + ".fasta")
    return sequence_set


def extract_based_country(fastafilename, country, path): 
    """extract all sequences for a specific country"""
    # if not os.path.exists(path):
    #     os.makedirs(path)
    handle = open(path + "sequences_" + country + ".fasta","w")
    for record in SeqIO.parse(fastafilename, "fasta"):
        try:
            str(record.id + record.description).index(country)
            print(record.id + record.description)
            SeqIO.write(record, handle,"fasta")
        except:
            continue
    handle.close()
        
        
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
    

def get_records_id_fasta(fastafilename): 
    index_set = []
    for record in SeqIO.parse(fastafilename, "fasta"):
        index_set.append(record.id)
    pd.DataFrame(index_set).to_csv("/alina-data0/Fatemeh/msa_0511_strain.csv")
    return index_set

    
    
def replace_space_in_ids(fasta_file):
    """replace the spaces in the sequences id with _"""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    handle = open(str(fasta_file).split('.fasta')[0] + "_space_fixed.fasta","w")
    for record in records:
        my_seq = SeqRecord(record.seq, id = record.description)
        SeqIO.write(my_seq, handle, "fasta") 
    handle.close()
    
    
def get_virus_names_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    ids = []
    for record in records:
        id_ = record.id + record.description
        ids.append("hCoV-19/" + id_.split("hCoV-19/")[-1].split('|')[0]) 
    return pd.DataFrame(ids, columns = ['Virus name'])

    
def get_epi_id_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    ids = []
    for record in records:
        id_ = record.id + record.description
        ids.append("EPI" + id_.split("EPI")[-1].split('|')[0]) 
    return pd.DataFrame(ids, columns = ['Accession ID'])


def seperate_metadata_months_continent(metadata):
    tsv_read = pd.read_csv(metadata, sep='\t')
    tsv_read = tsv_read[['Accession ID', 'Collection date', 'Virus name', 'Location', 'Host']]
    tsv_read = tsv_read.loc[tsv_read['Host'] == 'Human']
    location_list = []
    for location in tsv_read["Location"]:
        continent = location.split('/')[0].strip()
        location_list.append(continent) # + '_'+ country)
    tsv_read['continent'] = location_list
    dates = tsv_read['Collection date']
    tsv_read['Collection date m_y'] = pd.to_datetime(dates).dt.to_period('m')
    continent = ['Asia'] #North America', 'South America', 'Europe', 'Oceania', 'Africa']
    months =  tsv_read['Collection date m_y'].unique()
    for location in continent:
        df = tsv_read.loc[tsv_read['continent'] == location]
        for month in months:
            df.loc[df['Collection date m_y'] == month].to_csv('metadata_m_y/' + \
                                                              location + '_' + str(month) + '.csv')


if __name__ == '__main__':
    """continent metadata and path to fasta files for countries of that continent"""
    fasta_file = sys.argv[1] 
    # location = fasta_file.split('.fasta')[0].split('_')[-1]
    # path = '/alina-data0/Fatemeh/Epistasis/MSA/msa_1126/continents_metadata/' + location + '/'
    # dates = pd.date_range('2019-11-10','2021-11-26', freq='MS').strftime("%Y-%m").tolist()
    # for date in dates:
    extract_based_country(fasta_file, 'EPI_ISL_3808848', fasta_file)#'2021-04')
    # metadata = '/Users/fatemehmohebbi/Downloads/metadata_tsv_2021_12_09/metadata_2021_212_09.tsv'
    # seperate_metadata_months_continent(metadata)
    # extract_based_str_virus_name(df,file_path, file_path.split('.fasta')[0] + 'canada')
    # list_ = extract_based_submission_date(pd.to_datetime('12/1/2020'), metadata) #mm/dd/yyyy
    # list_[0].to_csv(str(metadata).split(".csv")[0] + '12_1_2020.csv')
    
