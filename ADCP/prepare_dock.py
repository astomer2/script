import os 
import sys
sys.path.append(os.path.abspath('.'))
from utils_peptide.pdb_fixer import fixed_pdb_file


def prepare_structure(protein, peptide, rec, lig):
    """
    Prepare protein and peptide structures for docking.

    Args:
        protein (str): The path to the protein PDB file.
        peptide (str): The path to the peptide PDB file.
        rec (bool): Whether to prepare the protein structure for docking.
        lig (bool): Whether to prepare the peptide structure for docking.
    """
    pdb_file_path, pdb_file = os.path.split(protein)
    fixed_pdb = fixed_pdb_file( protein, pdb_file_path)
    os.chdir(pdb_file_path)
    if rec :
        #os.system("reduce " + protein + " > " + protein.split(".")[0] + "-H.pdb")
        os.system("prepare_receptor -r " + fixed_pdb.split(".")[0] + ".pdb"+ " -o " + fixed_pdb.split(".")[0] + ".pdbqt")
    elif lig :
        #os.system("reduce " + peptide + " > " + peptide.split(".")[0] + "-H.pdb")
        os.system("prepare_ligand -l " + fixed_pdb.split(".")[0] + "-H.pdb" + " -o " + fixed_pdb.split(".")[0] + ".pdbqt")

    os.system("agfr -r " + fixed_pdb.split(".")[0] + ".pdbqt -b receptor -s 1 -p all -asv 1.1 -ps -o "+ fixed_pdb.split("/")[-1].split(".")[0])       
     


if __name__ == "__main__":
    protein = "/mnt/nas1/lanwei-125/TGFbR2/dock_prepare/TBR2.pdb"
    peptide = ""
    rec = True
    lig = False
    prepare_structure(protein, peptide, rec, lig)