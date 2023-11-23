
import os

protein_file = "/mnt/nas1/lanwei-125/IL8/v4/IL8dimer.txt"
peptides_file = "/mnt/nas1/lanwei-125/IL8/v4/seq.txt"
output_dir = "/mnt/nas1/lanwei-125/IL8/v4/AF/input/"

# 读取蛋白序列
with open(protein_file) as f:
    protein_seq = ''
    for line in f:
        if not line.startswith('>'):
            protein_seq += line.strip()

# 读取多肽序列  
with open(peptides_file) as f:
    peptide_seqs = []
    for line in f:
        if not line.startswith('>'):
            peptide_seqs.append(line.strip())

# 为每个多肽序列生成新的fasta文件
for peptide in peptide_seqs:
    outfile = os.path.join(output_dir, peptide + '_dimer.fasta')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(outfile, 'w') as f:
        f.write('>'+ peptide+'\n')
        f.write(protein_seq + ':\n')
        f.write( peptide)


print('Done. Output files written to', output_dir)


