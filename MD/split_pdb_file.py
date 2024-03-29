import os 
from genericpath import exists
from pathlib import Path
import shutil
import sys
sys.path.append(os.path.abspath('.'))
from utils_peptide.pdb_fixer import fixed_pdb_file

def make_dirs(file_path):
    for files in os.listdir(file_path):
        if files.endswith('.pdb'):
            filename = files.split(".")[0]
            source_path = os.path.join(file_path, files)
            dest_path = os.path.join(file_path, filename)
            os.umask(0)
            os.makedirs(dest_path, exist_ok = True, mode=0o777)
            shutil.copy(source_path, dest_path)

def read_pdb_chain(workdir):
    chains = {}
    workdir = Path(workdir)
    for pdb_name in [f for f in workdir.glob("*.pdb")]:
        if pdb_name.stem != 'peptide' and pdb_name.stem != 'protein':
            input_pdb = os.path.join(workdir, pdb_name)

            with open(input_pdb) as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain_id = line[21]
                        chains.setdefault(chain_id, []).append(line)
    return chains



def write_pdb_chain(chains, peptide, protein):
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

def split_pdb_run(workpath):
    make_dirs(workpath)
    for file in Path(workpath).iterdir():
        if file.is_dir():
            chains = read_pdb_chain(file)
            if len(chains) != 0:
                protein = os.path.join(file, 'protein.pdb')
                peptide = os.path.join(file, 'peptide.pdb')
                write_pdb_chain(chains, peptide, protein)
                fixed_pdb_file(peptide)
                fixed_pdb_file(protein)

if __name__ == '__main__':
    workpath = '/mnt/nas1/lanwei-125/MC5R/MD/AFD/'
    split_pdb_run(workpath)


'''

def main(workdir):
    protein = os.path.join(workdir, 'protein.pdb')
    peptide = os.path.join(workdir, 'peptide.pdb')
    chains = read_pdb_chain()
    write_pdb_chain(chains, peptide, protein)

if __name__ == '__main__':
    workdir = '/mnt/nas1/lanwei-125/IL8/v4/MD/HPSHFHG_monomer'
    main(workdir)

'''