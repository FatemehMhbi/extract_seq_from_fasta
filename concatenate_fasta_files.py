# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:36:37 2021

@author: Fatemeh
"""

import glob, os,sys
from Bio import SeqIO
import pandas as pd
#from utils import extract_strain_name, extract_EPI_id, extract_collection_date


def concate_fasta_files(fasta_files_path, pref_format):
    records = []
    for file in glob.glob(os.path.join(fasta_files_path, pref_format)):
        records = records + list(SeqIO.parse(file, "fasta"))
    handle = open(fasta_files_path + "all_concatenated.fasta", "w")
    SeqIO.write(records, handle, "fasta")
    handle.close()
    
    
def get_meta_from_fasta(seq_file, ref_file, id_name):
    sequences = list(SeqIO.parse(seq_file, "fasta"))
    ids = []
    for seq in sequences:
        ids.append(seq.id + seq.description)
    collection_date = extract_collection_date(ids)
    if id_name == 'Virus name':
        ids = extract_strain_name(ids)
    elif id_name == 'Accession ID':
        ids = extract_EPI_id(ids)
    metadata = pd.DataFrame({id_name: ids, 'Collection date': collection_date})
    
    
    
if __name__ == '__main__':
    path = sys.argv[1]
    concate_fasta_files(path, '*samples.fasta')