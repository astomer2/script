import csv
import pandas as pd
from pandas import DataFrame
from pathlib import Path
from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO
import sys 
import os
sys.path.append(os.path.abspath("."))

done=Path("/mnt/nas1/lanwei-125/TGFbR2/GA_generator/done.csv")
to_be_do=Path('/mnt/nas1/lanwei-125/TGFbR2/GA_generator/to_be_done.csv')
v4 = Path('/mnt/nas1/lanwei-125/TGFbR2/GA_generator/v4/cpp-classifier_prediction/cpp-generated_TGFbR2_prediction_all_pos.csv')
v5=Path('/mnt/nas1/lanwei-125/TGFbR2/GA_generator/v5/cpp-classifier_prediction/cpp-generated_TGFbR2_prediction_all_pos.csv')
df = pd.read_csv(done)
df1=pd.read_csv(to_be_do,header=None )
df2 = pd.read_csv(v4)
df2=df2['Sequence']
df3 = pd.read_csv(v5)
df3=df3['Sequence']

A =df.values.tolist()
B =df1.values.tolist()
C=df2.values.tolist()
D = []
for i in df3.values.tolist():
    if len(i) >5:
        D.append(i)
        print(i,len(i))

A =[a for b in A for a in b]
B =[a for b in B for a in b]
A=set(A)
B=set(B)
C=set(C)
D=set(D)
with open(to_be_do,"a+") as f:
    writer = csv.writer(f)
    for i in D-C-B-A:
        writer.writerow([i])