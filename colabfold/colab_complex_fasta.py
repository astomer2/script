
import os



def read_protein_fasta(file):
    with open(file) as f:
        protein_seq = ''
        for line in f:
            if not line.startswith('>'):
                protein_seq += line.strip()
    return protein_seq

def read_peptide_fasta(file):
    with open(file) as f:
        peptide_seqs = []
        for line in f:
            if not line.startswith('>'):
                peptide_seqs.append(line.strip())
    return peptide_seqs

# 为每个多肽序列生成新的fasta文件
def write_peptide_fasta(peptide_seqs, protein_seq, output_dir):
    for peptide in peptide_seqs:
        outfile = os.path.join(output_dir, peptide + '.fasta')
        if not os.path.exists(output_dir):
            os.umask(0)
            os.makedirs(output_dir)
        with open(outfile, 'w') as f:
            f.write('>'+ peptide+'\n')
            f.write(protein_seq + ':\n')
            f.write( peptide)

def main(protein_file, peptides_file, output_dir):
    protein_seq = read_protein_fasta(protein_file)
    peptide_seqs = read_peptide_fasta(peptides_file)
    write_peptide_fasta(peptide_seqs, protein_seq, output_dir)
    print('Done. Output files written to', output_dir)

if __name__ == '__main__':
    protein_file = "/mnt/sdc/lanwei/TGF/CF_relax/input/TBR2.fasta"
    peptides_file = "/mnt/sdc/lanwei/TGF/CF_relax/input/seq.txt"
    output_dir = "/mnt/sdc/lanwei/TGF/CF_relax/input/"
    main(protein_file, peptides_file, output_dir)

