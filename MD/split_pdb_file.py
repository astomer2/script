import os 

def read_pdb_chain():
    chains = {}
    for pdb_name in [f for f in os.listdir(workdir) if f.endswith('.pdb')]:
        input_pdb = os.path.join(workdir, pdb_name)

        with open(input_pdb) as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    chains.setdefault(chain_id, []).append(line)
    return chains



def write_pdb_chain(chains):
    chain_ids = len(chains)
    if chain_ids >= 2:
        all_chains = list(chains.values())

        # 设置最后一个值为peptide，其余为protein
        peptide_lines = all_chains[-1]
        protein_lines = all_chains[:-1]

        with open(peptide, 'w') as f:
            f.write(''.join(peptide_lines))

        with open(protein, 'w') as f:
            for lines in protein_lines:
                f.write(''.join(lines))
                f.write('TER\n')  # 添加TER行，表示链的结束

def main():

    chains = read_pdb_chain()
    write_pdb_chain(chains)
'''
if __name__ == '__main__':
    workpath = '/mnt/nas1/lanwei-125/IL8/v4/MD/'
    
    for file_name in os.listdir(workpath):
        workdir = os.path.join(workpath, file_name)
        if os.path.isdir(workdir):
            protein = os.path.join(workdir, 'protein.pdb')
            peptide = os.path.join(workdir, 'peptide.pdb')
            main()

'''

if __name__ == '__main__':
    workdir = '/mnt/nas1/lanwei-125/IL8/v4/MD/HPSHFHG_monomer'
    protein = os.path.join(workdir, 'protein.pdb')
    peptide = os.path.join(workdir, 'peptide.pdb')
    main()