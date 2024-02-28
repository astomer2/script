import os
from pprint import pprint
from pathlib import Path
import pandas as pd
import sys 
sys.path.append(os.path.abspath("."))
from utils_peptide.peptide_properties import calc_physico_chemical_properties
from utils_peptide.properties_from_structure.pep_structural_properties import calc_structural_properties

hPEP= Path('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence_best_pose/cluster/cluster_2')

A= []
for pdb in hPEP.glob('*.pdb'):
    A.append(pdb.stem)
len(A)
B= [seq for seq in A if seq.count("C")==2]
len(B)

df1 = pd.read_csv('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence/results.csv')
df2=df1[df1['sequence'].isin(B)]
df2.rename(columns={'sequence':'Sequence'},inplace=True)
df2 = df2.reset_index(drop=True)
df3 = calc_physico_chemical_properties(df2)

import shutil
os.umask(0)

os.makedirs('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence_best_pose/cluster/filter_structure/', exist_ok=True)
for seq in df3['Sequence']:
    shutil.copy(hPEP / (seq + '.pdb'), '/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence_best_pose/cluster/filter_structure/')
struct=Path('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence_best_pose/cluster/filter_complex/')
df4=calc_structural_properties(struct)
df5=pd.concat([df3,df4],axis=1)
df5.to_csv('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep-hpepdock_filter_structure.csv',index=False)