import os 
from genericpath import exists
from pathlib import Path
import shutil
import sys
import argparse
sys.path.append(os.path.abspath('.'))
from utils_peptide.pdb_fixer import fixed_pdb_file

def make_dirs(dir_path: str or Path): 
    dir_path = Path(dir_path)
    for files in dir_path.glob("*.pdb"):
        if files.is_file():
            new_dir= files.with_suffix("")
            os.umask(0)
            os.makedirs(files.with_suffix(""), exist_ok = True, mode=0o777)
            shutil.copy(files, new_dir)

def read_pdb_chain(workdir) -> dict:
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
    if len(chains) >= 2:
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
    return 

def split_pdb_run(workpath) -> None:
    workpath = Path(workpath)
    make_dirs(workpath)
    for dir in Path(workpath).iterdir():
        if dir.is_dir():
            chains = read_pdb_chain(dir)
            if len(chains) != 0:
                protein = dir/ 'protein.pdb'
                peptide = dir/ 'peptide.pdb'
                write_pdb_chain(chains, peptide, protein)
                fixed_pdb_file(peptide)
                fixed_pdb_file(protein)



if __name__ == '__main__':
    workdir = '/mnt/nas1/lanwei-125/PRLR/MD/pos1/'
    split_pdb_run(workdir)


'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='split pdb file')
    parser.add_argument('-i','--workpath', type=str, help='the path of workdir')
    args = parser.parse_args()
    workpath = args.workpath
    split_pdb_run(workpath)
'''
