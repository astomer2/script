import os 
import sys
from pathlib import Path
import shutil
import argparse
import time
sys.path.append(os.path.abspath('.'))
from utils_peptide.pdb_fixer import fixed_pdb_file


def prepare_structure(protein: str or Path, peptide: str or Path, trg_file: Path) -> None:  
    """
    Prepare the structure for molecular docking by copying the protein and peptide files to a new directory, 
    fixing the protein and peptide PDB files, and generating the receptor and ligand files in the PDBQT format.
    
    Parameters:
    - protein (str or Path): The path to the protein file.
    - peptide (str or Path): The path to the peptide file.
    
    Returns:
    - trg_file (Path): The path to the .trg file generated during the structure preparation.
    """

    protein = Path(protein)
    os.makedirs(protein.with_suffix(""), exist_ok=True)
    shutil.copy(protein, protein.with_suffix(""))
    pdb_file_path = protein.with_suffix("")
    os.chdir(pdb_file_path)
    fixed_protein_pdb = fixed_pdb_file(pdb_file_path/protein.name)
    if protein is not None:
        os.system("prepare_receptor -r " + fixed_protein_pdb.name+ " -o " + fixed_protein_pdb.stem + ".pdbqt") 
        if peptide is not None:
            peptide = Path(peptide)
            shutil.copy(peptide, protein.with_suffix(""))
            fixed_peptide_pdb = fixed_pdb_file(pdb_file_path/peptide.name)
            os.system("prepare_ligand -l " + fixed_peptide_pdb.name + " -o " + fixed_peptide_pdb.stem + ".pdbqt")
            os.system("agfr -r " + fixed_protein_pdb.stem + ".pdbqt -l " + fixed_peptide_pdb.stem + ".pdbqt -s 1 -p all -asv 1.1 -ps -o "+ fixed_protein_pdb.stem)
        else:
            os.system("agfr -r " + fixed_protein_pdb.stem + ".pdbqt -b receptor -s 1 -p all -asv 1.1 -ps -o "+ fixed_protein_pdb.stem)
        shutil.copy(pdb_file_path/(fixed_protein_pdb.stem+'.trg'), trg_file)
        #shutil.rmtree(pdb_file_path)
    return 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="prepare config files for docking")
    parser.add_argument("-p", "--protein", type=str, help="protein file path")
    parser.add_argument("-l", "--peptide", type=str, default=None, help="peptide file path")
    parser.add_argument("-o", "--trg_file", type=str, help="output .trg file path")

    args = parser.parse_args()
    protein = args.protein
    peptide = args.peptide
    trg_file = args.trg_file
    prepare_structure(protein, peptide, trg_file)