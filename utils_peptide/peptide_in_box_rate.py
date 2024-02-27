import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import PDB

sys.path.append(os.path.abspath("."))

from utils_comm.log_util import logger
from utils_peptide.pdb_utils import get_sequence_from_chain


def get_ca_coordinates(structure_path, chain_id, residue_number):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', structure_path)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == residue_number and 'CA' in residue:
                        return residue['CA'].get_coord()
    return None

def calculate_percentage_inside_box(box, coordinates):
    inside_box = np.all((coordinates >= box[0]) & (coordinates <= box[1]), axis=1)
    return np.sum(inside_box) / len(coordinates)

def calc_in_box_rate(pdb_file_path, target_chain_id, target_residue_number,other_chain_id):
    ca_coords = get_ca_coordinates(pdb_file_path, target_chain_id, target_residue_number)
    if ca_coords is not None:
        box = np.array([ca_coords - 20, ca_coords + 20])
        other_chain_id = "B"
        other_chain_coordinates = []
        structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb_file_path)
        for model in structure:
            for chain in model:
                if chain.id == other_chain_id:
                    for residue in chain:
                        if 'CA' in residue:
                            other_chain_coordinates.append(residue['CA'].get_coord())
        percentage_inside_box_rate = calculate_percentage_inside_box(box, np.array(other_chain_coordinates))
        tartget_chain_seq = get_sequence_from_chain(pdb_file_path, other_chain_id)
        return tartget_chain_seq, percentage_inside_box_rate


def multiple_calc_in_box_rate(
    pdb_dir_path, target_chain_id, target_residue_number, other_chain_id, saved=True
):
    result = []
    pdb_dir_path = Path(pdb_dir_path)
    for pdb_file_path in pdb_dir_path.glob('*.pdb'):
        tartget_chain_seq, percentage_inside_box_rate = calc_in_box_rate(
            str(pdb_file_path), target_chain_id, target_residue_number, other_chain_id
        )
        append_result = [pdb_file_path.stem, tartget_chain_seq, percentage_inside_box_rate]
        result.append(append_result)
    df = pd.DataFrame(result, columns=['file_name', 'target_chain_seq', 'percentage_inside_box_rate'])
    if saved:
        df.to_csv(pdb_dir_path/'in_box_rate.csv', index=False)
        return df
    else:
        return df


def keep_complexes_binded_at_correct_location(
    complex_pdb_dir,
    target_chain_id,
    target_residue_number,
    other_chain_id="B",
    threshold_ratio=0.6,
):
    complex_pdb_dir = Path(complex_pdb_dir)
    all_df = multiple_calc_in_box_rate(
        complex_pdb_dir,
        target_chain_id=target_chain_id,
        target_residue_number=target_residue_number,
        other_chain_id=other_chain_id,
        saved=True,
    )
    correct_binded_pdb_dir = complex_pdb_dir / "correct_binded_pdbs"
    correct_binded_pdb_dir.mkdir(parents=True, exist_ok=True)
    result_file = correct_binded_pdb_dir / "correct_binded_pdbs.csv"
    df = all_df[all_df["percentage_inside_box_rate"] >= threshold_ratio]
    all_pdb_num = len(all_df)
    correct_binded_pdb_num = len(df)
    correct_binded_pdb_ratio = correct_binded_pdb_num / all_pdb_num
    logger.info(
        f"{all_pdb_num = }, {correct_binded_pdb_num = }, {correct_binded_pdb_ratio = :.2f}"
    )
    df.to_csv(result_file, index=False)
    for file_stem in df["file_name"]:
        pdb_filename = f"{file_stem}.pdb"
        shutil.copy(
            complex_pdb_dir / pdb_filename, correct_binded_pdb_dir / pdb_filename
        )
        json_file = complex_pdb_dir / f"{file_stem}.json"
        if json_file.is_file():
            shutil.copy(json_file, correct_binded_pdb_dir / json_file.name)


if __name__ == '__main__':
    pdb_dir_path = '/mnt/nas1/lanwei-125/MC5R/MD/AFD/'
    target_chain_id = 'A'  # 指定的受体链的ID
    target_residue_number = 89  # 指定的受体氨基酸残基号
    other_chain_id='B'
    multiple_calc_in_box_rate(pdb_dir_path, target_chain_id, target_residue_number,other_chain_id,saved=True)
