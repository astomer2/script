import os 
import sys
from pathlib import Path
import logging 
import time
import csv
import pandas as pd
from pandas import DataFrame
from pathlib import Path
from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO
sys.path.append(os.path.abspath("."))
from utils_peptide.pdb_utils import get_sequence_from_chain
from utils_peptide.peptide_utils import get_specify_aa_position

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
now_times = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def read_seq_pdb(pdb_path:str or Path) -> str:
    parse = PDBParser()
    pdb_path=Path(pdb_path)
    structure = parse.get_structure('X', pdb_path)
    resname_str = ' '.join([seq.get_resname() for model in structure for chain in model for seq in chain])
    resname_str_txt = pdb_path.with_suffix('.txt')
    chain_seq= get_sequence_from_chain(pdb_path)
    with open(resname_str_txt, 'w') as f:
        f.write(resname_str)
    return resname_str_txt,chain_seq

def code_distribute(np_nums):
    director = 1
    manager = np_nums/8
    worker = np_nums-director-manager
    return director,manager,worker

def calc_pnear(pdb_path, np_nums, cys=False, cyx=False):
    pdb_path=Path(pdb_path)
    logger.info(f'{now_times}: processing {pdb_path}')

    resname_str_txt,chain_seq = read_seq_pdb(pdb_path)
    cyclization_type = None
    A,B,C = code_distribute
    logger.info(f"{now_times}: calculate {chain_seq} using {np_nums} process")

    f=open(pdb_path.parent/'pnear.flags', 'w')
    f.write(f"-cyclic_peptide:MPI_processes_by_level {A} {B} {C}\n")
    f.write("-cyclic_peptide:MPI_batchsize_by_level 400 4 \n")
    f.write("-cyclic_peptide:MPI_output_fraction 1 \n")
    f.write("-nstruct 50000 \n")
    f.write(f"-cyclic_peptide:sequence_file {resname_str_txt} \n")
    if cys and chain_seq.count("C") >= 2:

        cyclization_type = 'terminal_disulfide'
        C_position=get_specify_aa_position(chain_seq, 'C')
        f.write(f"-cyclic_peptide:exclude_residues_from_rms {C_position[0]+1}  {C_position[1]+1} \n")
        logger.info(f"{now_times}: {pdb_path} has two cys, use terminal_disulfide")
    if cyx:
        cyclization_type = 'n_to_c_amide_bond'

    f.write(f"-cyclic_peptide:cyclization_type {cyclization_type} \n")
    f.write("-cyclic_peptide:genkic_closure_attempts 1000 \n")
    f.write("-cyclic_peptide:min_genkic_hbonds 1 \n")
    f.write("-cyclic_peptide:MPI_sort_by energy \n")
    f.write("-cyclic_peptide:compute_pnear_to_this_fract 1 \n")
    f.write("-cyclic_peptide:fast_relax_rounds 20 \n")
    f.write("-mute all \n")
    f.write("-ex1 \n")
    f.write("-ex2 \n")
    f.write("-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary \n")
    f.write(f"-in:file:native {pdb_path} \n")
    f.write(f"-out:file:silent {pdb_path.with_suffix('.out')} \n")
    f.close()
    log_file= pdb_path.with_suffix('.log')
    os.system(f"mpirun.openmpi --use-hwthread-cpus -np {np_nums} simple_cycpep_predict.mpi.linuxgccrelease @{pdb_path.parent/'pnear.flags'}> {log_file} 2>&1")

    with open(log_file, 'r') as f:
        lines = f.readlines()
        found_first = False  
        for line in lines:
            if line.startswith("MPI_worker_node"):
                if found_first:  
                    next_line = (lines[lines.index(line)+1])
                    best_pnear=next_line.split()[2]
                    print(best_pnear)
                    break  
                found_first = True  
    logger.info(f"{now_times}:{chain_seq} best_pnear is {best_pnear}")
    with open(pdb_path.with_suffix('.best_pnear'), 'w') as f:
        f.write(chain_seq,best_pnear)


if __name__ == "__main__" :
    pdb_path="/mnt/nas1/lanwei-125/test/pnear/CAEVWIEHC.pdb"
    np_nums=20
    cys= True
    calc_pnear(pdb_path ,np_nums, cys)
