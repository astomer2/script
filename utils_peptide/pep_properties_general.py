import os
import random
import sys
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd
from pandas import DataFrame

sys.path.append(os.path.abspath("."))
from utils_comm.log_util import logger
from utils_peptide.peptide_properties import calc_physico_chemical_properties
from utils_peptide.properties_from_structure.pep_structural_properties import (
    calc_structural_properties,
)


SEED = 0
random.seed(SEED)


def calc_pep_properties(
    input_sequences: Optional[Union[DataFrame, List[str], str, Path]] = None,
    input_pdb_path: Optional[Union[str, Path]] = None,
    output_file=None,
) -> DataFrame:
    """Calculate peptide properties from sequence and structure

    Args:
        input_sequences: sequences of peptides
        input_pdb_path: Path to the single input PDB file or directory containing multiple PDB files
        output_file: output file

    If input_sequences and input_pdb_files_dir are both not None, please assert the seq
    number and input_pdb_files number are the same, and their orders are also the same!
    If possible, uses the sequence as the pdb file name.

    Returns:
        DataFrame: properties of peptides
    """
    if input_sequences is not None:
        df = calc_physico_chemical_properties(input_sequences)
    else:
        df = None
    if input_pdb_path is not None:
        # validate the number of sequences and pdb files.
        if df is not None:
            input_pdb_path = Path(input_pdb_path)
            if input_pdb_path.is_file() and input_pdb_path.suffix == ".pdb":
                assert (
                    len(df) == 1
                ), f"The number of sequences is {len(df)}, but there is 1 PDB file, not equal!"
            elif input_pdb_path.is_dir():
                pdb_file_num = len(list(input_pdb_path.glob("*.pdb")))
                assert (
                    len(df) == pdb_file_num
                ), f"The number of sequences {len(df)} is not equal to the number of PDB files {pdb_file_num}!"

        df_structural = calc_structural_properties(input_pdb_path)
        if df is not None:
            assert len(df) == len(df_structural)
            df = pd.concat([df, df_structural], axis=1)
        else:
            df = df_structural
    if input_sequences is None and input_pdb_path is None:
        raise ValueError("Both input_sequences and input_pdb_files are None")
    if output_file is not None and df is not None:
        df.to_csv(output_file, index=False)
    if df is not None:
        logger.info("df.columns %s, df.head()\n%s", df.columns.tolist(), df.head())
    return df  # type: ignore


if __name__ == "__main__":
    calc_pep_properties(
        input_sequences=["AFDGHLKI", "DLLRTK", 'AFDGHL', 'KKLQRS', 'AFDLKI', 'KKLLRTK'],
        input_pdb_path="/mnt/nas1/lanwei-125/test/",
    )
    logger.info("end")
