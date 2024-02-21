import os 
import pandas as pd
from pathlib import Path
import csv
df1=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/MRAP2/motif_4/MRAP2_4_filter.csv')
df2=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/POMC/motif_4/POMC_4_filter.csv')
df3=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/MRAP/motif_4/MRAP_4_filter.csv')
df4=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/ASIP/motif_4/ASIP_4_filter.csv')
df5=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/AGRP/motif_4/AGRP_4_filter.csv')

df6=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/MRAP2/motif_5/motif_5_5_filter.csv')
df7=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/POMC/motif_5/motif_5_5_filter.csv')
df8=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/MRAP/motif_5/motif_5_5_filter.csv')
df9=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/ASIP/motif_5/motif_5_5_filter.csv')
df10=pd.read_csv('/mnt/nas1/lanwei-125/MC5R/cPEP/AGRP/motif_5/motif_5_6_filter.csv')
all_sequences = df1['Sequence'].tolist() + df2['Sequence'].tolist() + df3['Sequence'].tolist() + df4['Sequence'].tolist() + df5['Sequence'].tolist() + df6['Sequence'].tolist()+ df7['Sequence'].tolist()+ df8['Sequence'].tolist()+ df9['Sequence'].tolist()+ df10['Sequence'].tolist()

 
all_sequences=list(set(all_sequences))
with open('/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence//Sequence.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    writer.writerow(['Sequence'])
    for seq in all_sequences:
        writer.writerow([seq])