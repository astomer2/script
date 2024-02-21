import os 

def read_receptor(protein_path, old_chain_id,new_chain_id):
    receptor = []
    with open(protein_path) as f:
        for line in f:
            #抓取ATOM行，用正则取出残基名称进行修改
            if line.startswith('ATOM'):
                residue = line[17:20]
                if residue == 'CYX':
                    line = line[:17] + 'CYS' + line[20:]
                
                chain_id = line[20]
                if chain_id != old_chain_id:
                    line = line[:21] + new_chain_id + line[22:]

                occupancy = line[54:60]
                if occupancy != '0.00':
                    line = line[:56] + '1.00' + line[60:]
                receptor.append(line)
    return receptor

def read_ligand(peptide_dir_path, complex_output_path):
    os.umask(0)
    os.makedirs(complex_output_path, exist_ok=True)
    for root, dirs, files in os.walk(peptide_dir_path):
        for file in files:
            if file.endswith(".pdb"):
                peptide_lines =[]
                with open(os.path.join(root, file)) as f:

                    lines = f.readlines()
                    for line in lines:
                        if line.startswith('ATOM'):
                            peptide_lines.append(line)

                complex_pdb = os.path.join(complex_output_path, file)
                with open(complex_pdb, 'w') as f:
                    f.writelines(new_pdb)
                    f.writelines('\nTER\n')
                    f.writelines(peptide_lines)
                    f.writelines("END")

if __name__ == '__main__':

    protein_path = '/mnt/nas1/lanwei-125/MC5R/dock/structure/MC5R.pdb'
    peptide_dir_path = '/mnt/nas1/lanwei-125/MC5R/dock/pos/'
    complex_output_path = '/mnt/nas1/lanwei-125/MC5R/dock/pos/'
    old_chain_id = "A"
    new_chain_id = "A"

    new_pdb = read_receptor(protein_path, old_chain_id, new_chain_id)
    read_ligand(peptide_dir_path,complex_output_path)