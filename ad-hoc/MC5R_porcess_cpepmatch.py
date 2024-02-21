import csv
import pandas as pd
from pandas import DataFrame
from pathlib import Path
from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO
import sys 
import os
sys.path.append(os.path.abspath("."))


three_to_one_mapping = {
    'CYS': 'C',
    'HIS': 'H',
    'PRO': 'P',
    'GLN': 'Q',
    'VAL': 'V',
    'ASP': 'D',
    'SER': 'S',
    'PHE': 'F',
    'GLY': 'G',
    'ARG': 'R',
    'GLU': 'E',
    'TYR': 'Y',
    'ILE': 'I',
    'LEU': 'L',
    'ASN': 'N',
    'MET': 'M',
    'TRP': 'W',
    'ALA': 'A',
    'THR': 'T',
    'LYS': 'K',
    "ACE": "" ,
    "NH2": ""
}

def extract_three_letter_code(peptide_structure):
    sequences = []

    with open(peptide_structure, "r", encoding="utf-8") as f:
        lines = f.readlines()
    sequence = []
    residue_ids = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM") :
            res_id = line[22:27].strip()
            if res_id not in residue_ids:
                residue_ids.append(res_id)
                res_name = line[17:20]
                if res_name not in ["HOH", "NH2", "ACE", "NME"]:
                    sequence.append(res_name)
    sequence_str = " ".join(sequence)
    sequences.append(sequence_str)
    sequence_lengths = len(sequence)
    return sequences, sequence_lengths


CPEP=Path('/mnt/nas1/lanwei-125/MC5R/cPEP/ASIP/motif_6/')

cpep=[]
for i in CPEP.glob('match*.pdb'):
    if "NotMutated" not in i.name:
        sequences, sequence_lengths = extract_three_letter_code(i)
        
        if 5 < sequence_lengths < 11:
            converted_sequences_list = []  # 创建一个空列表来存储转换后的序列
            for sequence in sequences:
                converted_sequence = ''.join([three_to_one_mapping[code] for code in sequence.split() if code in three_to_one_mapping])
                converted_sequences_list.append(converted_sequence)
            
            # 使用 join 将列表中的所有序列连接在一起
            converted_sequences = ' '.join(converted_sequences_list)
            a = converted_sequences,len(converted_sequences),i.name
            cpep.append(a)
            #print(converted_sequences)

df = DataFrame(cpep)
df_unique = df.drop_duplicates(subset=0,keep='first')
df_unique.columns = ['Sequence','length','pdb_name']
df_unique.to_csv(CPEP/(str(CPEP.name)+"_6_filter.csv"),index=False)