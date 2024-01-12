import os
import sys

import shutil
import pandas as pd
import argparse
import logging
import time
from math import log, log10
from pandas import DataFrame
from pathlib import Path


'''
#按照pdb文件名称新建文件夹，注意需要是蛋白-肽复合体
def read_pdb():
    pdb_names = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith('.pdb'):
            pdb_names.append(file_name.split('.')[0])
            new_dir = os.path.join(dir_path, file_name.split('.')[0])
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            pdb_file = os.path.join(dir_path, file_name)
            shutil.copy(pdb_file, new_dir)
    return pdb_names
'''    
# 此脚本将与聚类结果偶联，如果没做或者没有聚类文件夹,则启用上面被注释的代码，记得修改下面的路径命名

#读取pdb文件所在的工作路径

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
now_times = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def read_pdb(dir_path: str):
    """
    Reads the provided directory path and returns a list of cluster directory names.

    Parameters:
        dir_path (str): The path to the directory containing cluster directories.

    Returns:
        List[str]: A list of cluster directory names.
    """
    cluster_dir_names = [
        dir_name for dir_name in Path(dir_path).iterdir()
        if dir_name.is_dir() and dir_name.name.split("_")[0] == "cluster"
    ]

    return [str(cluster_dir_name) for cluster_dir_name in cluster_dir_names]

#结构预处理，将残基命名，原子参数等设置正确，为配体添加受体，受体添加在配体之前。rosetta只能读取标注氨基酸命名，
#如果必须包含非标准氨基酸，则需要在flags文件中写明
def preprocess_pdb(dir_path: str, protein_pdb: str, cluster_dir_names: list):
    """
    Preprocesses a protein PDB file by extracting ATOM and TER lines.
    
    Args:
        dir_path (str): The directory path where the protein PDB file and cluster directories are located.
        protein_pdb (str): The path to the protein PDB file.
        cluster_dir_names (List[str]): A list of cluster directory names.
    
    Returns:
        None
    
    Raises:
        FileNotFoundError: If the protein PDB file or any of the cluster directories or files are not found.
    """
    protein_lines = []
    protein_chain = set()
    with open(protein_pdb) as f:
        protein_lines = [line for line in f if line.startswith("ATOM") or line.startswith("TER")]
        for line in protein_lines:
            if line.startswith('ATOM') :
                chains = line[21]
                protein_chain.add(chains)
    protein_chain_str = ','.join(protein_chain)

    peptide_chain = set()
    peptide_chain_str = ''
    for cluster_dir_name in cluster_dir_names:
        pdb_path = Path(dir_path) / cluster_dir_name

        for pdb_file in pdb_path.glob('*.pdb'):
            parts = pdb_file.stem.split("-")
            if len(parts) >= 2 and parts[1] == "ref":
                continue

            if pdb_file.suffix == '.pdb' and '-ref' not in pdb_file.stem:
                # Read peptide lines
                peptide_lines = []
                with open(pdb_file) as f:
                    peptide_lines = [line for line in f if line.startswith("ATOM")]
                    for line in peptide_lines:
                        if line.startswith('ATOM') :
                            chains = line[21]
                            peptide_chain.add(chains)
                peptide_chain_str = ','.join(peptide_chain)

                # Write to new pdb file
                os.makedirs(pdb_path / pdb_file.stem, exist_ok=True)
                new_pdb_path = pdb_path / pdb_file.stem / (pdb_file.stem + '-ref.pdb')
                with open(new_pdb_path, 'w') as f:
                    f.write(''.join(protein_lines))
                    if protein_lines[-1].startswith('END'):
                        #删除END语句
                        del protein_lines[-1]
                    if not protein_lines[-1].startswith('TER'):
                        f.write("\nTER\n")
                    f.write(''.join(peptide_lines))
                    f.write("\nTER\n")
                    f.write("END\n")

                logger.info(f'{now_times}: processing {pdb_file} to {new_pdb_path}')

    return protein_chain_str,peptide_chain_str

# 使用rosetta_flex_pepdock进行结构细化处理
# 首先预处理pdb文件，然后进行局部细化，输出打分文件refinement.score.sc
def make_flag(dir_path: str, cluster_dir_names: list, protein_chain: str, peptide_chain: str, np_nums: str):
    """
    Generates flags and runs FlexPepDocking for each cluster directory in the specified directory path.

    Args:
        dir_path (str): The path to the directory containing the cluster directories.
        cluster_dir_names (list): A list of cluster directory names.
        protein_chain (str): The protein chain.
        peptide_chain (str): The peptide chain.
        np_nums (int): The number of processes to run.

    Returns:
        None
    """
    np_nums = str(np_nums)
    protein_chain_id = str(protein_chain)
    peptide_chain_id = str(peptide_chain)
    for cluster_dir_name in cluster_dir_names:

        for file in (Path(dir_path) / cluster_dir_name).glob('*.pdb'):
            pdb_file = file.stem
            work_dir = Path(dir_path) / cluster_dir_name / pdb_file
            os.chdir(work_dir)

            for file in work_dir.glob('*-ref.pdb'):
                ref_pdb = file.stem
                logger.info(f'{now_times}: processing {file}')

                f=open('prepack.flags', 'w')
                f.write("-s "+ ref_pdb + ".pdb \n")
                f.write("-flexpep_prepack \n")
                f.write("-nstruct 1 \n")
                f.write("-ex1 \n")
                f.write("-ex2aro \n")
                f.write("-mute core.io.database \n")
                f.write("-mute core.util.prof \n")
                f.write("-mute protocols.jd2.JobDistributor \n")
                f.write("-scorefile pack.score.sc ")
                f.close()
                os.system("mpirun.openmpi -np "+ np_nums + " FlexPepDocking.mpi.linuxgccrelease @prepack.flags")

                f=open('refinement.flags', 'w')
                f.write("-s "+ ref_pdb + "_0001.pdb \n")
                #f.write("-native "+ peptide + ".pdb \n")
                f.write("-flexPepDocking:receptor_chain "+protein_chain_id+" \n")
                f.write("-flexPepDocking:peptide_chain "+peptide_chain_id+" \n")
                f.write("-pep_refine \n")
                f.write("-use_input_sc \n")
                f.write("-nstruct 250 \n")
                f.write("-scorefile refinement.score.sc \n")
                f.write("-ex1 \n")
                f.write("-ex2aro \n")
                f.write("-mute core.io.database \n")
                f.write("-mute core.util.prof \n")
                f.write("-mute protocols.jd2.JobDistributor ")
                f.close()
                os.system("mpirun.openmpi -np "+ np_nums + " FlexPepDocking.mpi.linuxgccrelease @refinement.flags")

# 从refinement.score.sc中选出最优结构
def take_candidate(cluster_dir_names: list, result_path: Path):
    """
    Takes a list of cluster directory names and a result path as input.

    Creates the result path directory if it doesn't exist.

    For each cluster directory name in the given list:
        - Iterates over all the pdb files in the cluster directory.
        - Parses the pdb file stem and constructs the work directory path.
        - Constructs the path to the score file in the work directory.
        - Reads the score file and extracts the headers and rows.
        - Creates a pandas DataFrame from the rows and headers.
        - Sorts the DataFrame by the 'reweighted_sc' column.
        - Retrieves the row with the minimum 'reweighted_sc' value.
        - Constructs the path to the candidate pdb file.
        - Copies the candidate pdb file to the result path.
        - Writes the candidate description and reweighted score to the result score file.

    Parameters:
    - cluster_dir_names (list): A list of cluster directory names.
    - result_path (str): The path to the result directory.

    Returns:
    - files path that contains the candidate pdb files.
    """
    result_path = Path(result_path)

    if not result_path.exists():
        result_path.mkdir()

    for cluster_dir_name in cluster_dir_names:
        for file in (Path(dir_path) / cluster_dir_name).glob('*.pdb'):
            pdb_file = file.stem
            work_dir = Path(dir_path) / cluster_dir_name / pdb_file  
            score_file = work_dir / 'refinement.score.sc'

            with open(score_file, 'r') as f:
                headers = None
                rows = []
                lines = f.readlines()
                for line in lines:
                    if line.startswith('SEQUENCE:'):
                        continue
                    elif line.startswith('SCORE:'):
                        if headers is None:
                            headers = line.split()[1:]
                        else:
                            rows.append(line.split()[1:])
                    else:
                        rows.append(line.split())
                df = pd.DataFrame(rows, columns=headers)
                df['reweighted_sc'] = df['reweighted_sc'].astype(float)
                df = df.sort_values(by='reweighted_sc', ignore_index=True)
                min_row = df.loc[df['reweighted_sc'].idxmin()]

            candidate_path = work_dir / f"{min_row['description']}.pdb"
            shutil.copy(candidate_path, result_path)

            result_score_file = result_path / 'candidate.txt'
            with open(result_score_file, 'a') as f:
                f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")
    return result_score_file

# 局部细化后处理
def process_refine_pdb(result_path):
    """
    Process the refine PDB files in the specified directory.

    This function iterates over all files in the given directory and filters out
    files with the '.pdb' extension. It then processes each of these files by
    extracting lines that start with 'ATOM', 'TER', or 'END' and writes them
    to a new file with the same stem as the original file.

    Parameters:
    - result_path (str): The path to the directory containing the refine PDB files.

    Returns:
    - None

    Raises:
    - FileNotFoundError: If the specified directory does not exist.
    """
    for file_entry in Path(result_path).iterdir():
        if file_entry.suffix == '.pdb':
            processed_pdb_path = Path(result_path) / file_entry.name
            candidate_pdb_path = Path(result_path) / f"{file_entry.stem.split('-')[0]}.pdb"
            candidate_pdb = []

            with open(processed_pdb_path) as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
                        candidate_pdb.append(line)

            with open(candidate_pdb_path, 'w') as f:
                f.write(''.join(candidate_pdb))

def extract_pdb(candidate_txt):
    
    candidate_txt = Path("")
    with open(candidate_txt, 'r') as f:
        lines = f.readlines()
        c = []
    for line in lines:
        a = line.split()[0].split("-")[0]
        b = line.split()[1]
        if log10(abs(float(b))) < 4:
            newline = a + " " + b 
            c.append(newline)
    df1 = pd.DataFrame(c)
    d = df1[0].str.split(' ', expand=True)
    df1['seq'] = d[0]
    df1['value'] = d[1]
    df1 = df1.drop(columns=[0])
    df1.sort_values(by=['value'], inplace=True,ascending=False)
    df1 = df1.reset_index(drop=True)
    df1 = df1.iloc[:,50:]
    for seq in df1['seq']:
        os.umask(0) 
        HPEP = candidate_txt.parent/(seq + ".pdb")
        hpep_md = "/mnt/nas1/lanwei-125/FGF5/disulfide/MD/HPEP/"
        os.makedirs(hpep_md, exist_ok=True)
        shutil.copy(HPEP, hpep_md)


def run_refinement(dir_path, result_path, protein_pdb, np_nums):
    cluster_dir_names = read_pdb(dir_path)
    protein_chain, peptide_chain = preprocess_pdb(dir_path, protein_pdb, cluster_dir_names)
    make_flag(dir_path, cluster_dir_names, protein_chain, peptide_chain, np_nums)
    result_score_file =take_candidate(cluster_dir_names, result_path)
    process_refine_pdb(result_path)
    #extract_pdb( result_score_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='after clustering, run refinement, use rosetta to refine the pdb, choose the best one as the final pdb')
    parser.add_argument("-i", "--dir_path", type=str, help='you cluster result in dir_path')
    parser.add_argument("-o", "--result_path", type=str, help='output refinement to result_path')
    parser.add_argument("-rec", "--protein_pdb", type=str, help='recepor protein_pdb')
    parser.add_argument("-n", "--np_nums", type=int, help='np_nums, max number of half cores processing')

    args = parser.parse_args()
    dir_path = args.dir_path
    result_path = args.result_path
    protein_pdb = args.protein_pdb
    np_nums = args.np_nums

    run_refinement(dir_path, result_path, protein_pdb, np_nums)