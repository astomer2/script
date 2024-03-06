import pandas as pd 
from pathlib import Path
import os
from Levenshtein import distance
from sklearn.cluster import OPTICS

ex=Path('/mnt/nas1/lanwei-125/test/新建 Microsoft Excel 工作表.xlsx')

df = pd.read_excel(ex,sheet_name='Sheet5')
seqs = df['peptide_seq'].to_list()

new_sequences = []

for sequence in seqs:
    if 'R' in sequence:
        r_index = sequence.index('R')
        new_sequence = sequence[r_index:] + sequence[:r_index]
        new_sequences.append(new_sequence)
    elif 'F' in sequence:
        f_index = sequence.index('F')
        new_sequence = sequence[f_index:] + sequence[:f_index]
        new_sequences.append(new_sequence)
    elif 'I' in sequence:
        i_index = sequence.index('I')
        new_sequence = sequence[i_index:] + sequence[:i_index]
        new_sequences.append(new_sequence)
    else:
        new_sequences.append(sequence)

df['new_peptide_seq'] = new_sequences

df.to_excel('/mnt/nas1/lanwei-125/test/新建 Microsoft Excel 工作表.xlsx', index=False)

length_dfs = {}

for length, group_df in df.groupby(df['new_peptide_seq'].str.len()):
    length_dfs[length] = group_df

for length, length_df in length_dfs.items():
    distances = []
    for i, seq1 in enumerate(length_df['new_peptide_seq']):
        row_distances = [distance(seq1, seq2) for seq2 in length_df['new_peptide_seq']]
        distances.append(row_distances)

    clustering =OPTICS(min_samples=2)
    clusters = clustering.fit_predict(distances)
    length_df['Cluster'] = clusters
    sorted_df = length_df.sort_values(by='Cluster', inplace=True)
    length_df.to_csv(ex.parent/f'length_{length}_all_clusters.csv', index=False)    