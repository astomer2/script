#TO DO ,Integration of three tools for calculating structural properties
import json
import math
import os
import random
import re
import shutil
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd
from pandas import DataFrame
from tqdm import tqdm
sys.path.append(os.path.abspath("."))
#from utils_peptide.peptide_properties import log_df_basic_info
from utils_comm.log_util import logger
from utils_peptide.properties_from_structure.caculate_potential_energy import calc_complex_energy,copy_and_create_directory,read_pdb_chain
from utils_peptide.properties_from_structure.calc_relative_SASA import calc_SASA
from utils_peptide.properties_from_structure.calc_residue_depth import interface_residues_depth

def structure_properties_to_dataframe(input_pdb_files: Union[str, Path], ) -> DataFrame:
    """
    Converts structure properties from PDB files to a pandas DataFrame.

    Args:
        input_pdb_files (Union[str, Path]): The path to the input PDB file.

    Returns:
        DataFrame: A DataFrame containing the structure properties.
    """
    input_pdb_files = Path(input_pdb_files)
    energy = calc_complex_energy(input_pdb_files)

    protein = input_pdb_files.with_suffix("")/"protein.pdb"
    peptide = input_pdb_files.with_suffix("")/"peptide.pdb"    

    protein_sasa = calc_SASA(protein)
    peptide_sasa = calc_SASA(peptide)
    complex_sasa = calc_SASA(input_pdb_files)
    relative_sasa = complex_sasa - (protein_sasa + peptide_sasa)
    #pdb_path, chain_id, res_id, RD[0], RD[1]
    #depth = interface_residues_depth(input_pdb_files)
    data = {
        'file_name': [energy.file_name],
        'raw_complex_energy': [energy.raw.complex_energy],
        'raw_protein_energy': [energy.raw.protein_energy],
        'raw_peptide_energy': [energy.raw.peptide_energy],
        'raw_diff_energy': [energy.raw.diff_energy],
        'minimized_complex_energy': [energy.minimized.complex_energy],
        'minimized_protein_energy': [energy.minimized.protein_energy],
        'minimized_peptide_energy': [energy.minimized.peptide_energy],
        'minimized_diff_energy': [energy.minimized.diff_energy],
        'protein_sasa': [protein_sasa],
        'peptide_sasa': [peptide_sasa],
        'complex_sasa': [complex_sasa],
        'relative_sasa': [relative_sasa]
    }

    df = pd.DataFrame(data)
    df.reset_index(drop=True, inplace=True)
    return df

def calc_structural_properties(
        input_pdb_path: Optional[Union[str, Path]] = None
) -> DataFrame :
    """
    Calculates the structural properties of protein-peptide complexes, NEED openmm

    Args:
        input_pdb_path (Optional[Union[str, Path]]): Path to the input PDB file or directory containing multiple PDB files. Defaults to None.

    Returns:
        DataFrame: DataFrame containing the calculated structural properties.

    Raises:
        ValueError: If the input PDB file type is not supported or if no valid protein-peptide complex PDB files are found in the given directory.
    """
    input_pdb_path = Path(input_pdb_path)
    if input_pdb_path.is_file():
        if input_pdb_path.suffix == ".pdb" and len(read_pdb_chain(input_pdb_path))>=2:
            df = structure_properties_to_dataframe(input_pdb_path)
            shutil.rmtree(input_pdb_path.with_suffix(""))
            return df
        else:
            raise ValueError(f"input_pdb_files type {input_pdb_path.suffix} not supported,should be a protein-peptide complex pdb file")
    elif input_pdb_path.is_dir():
        df_list = []
        for pdb_file in input_pdb_path.glob("*.pdb"):
            if pdb_file.is_file() and len(read_pdb_chain(pdb_file)) >= 2:
                logger.info(f"Processing {pdb_file}")
                df = structure_properties_to_dataframe(pdb_file)
                df_list.append(df)
                shutil.rmtree(pdb_file.with_suffix(""))
            else:
                raise ValueError(f"No valid protein-peptide complex pdb files found in {input_pdb_path}")
        integrated_df = pd.concat(df_list, axis=0, ignore_index=True)
        return integrated_df

    else:
        raise ValueError(f"input_pdb_files type {pdb_file.suffix} not supported, should be a folder containing protein-peptide complex pdb files")


if __name__ == "__main__":
    input_pdb_path = '/mnt/nas1/lanwei-125/MC5R/dock/pos/'
    df = calc_structural_properties(input_pdb_path)
    print(df)
    df.to_csv(f"{input_pdb_path}structure_properties.csv", index=False)
    print(f"Structural properties calculated and saved to {input_pdb_path}.csv")
