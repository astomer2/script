import os
from pdbfixer import PDBFixer
from pathlib import Path
from openmm.app import PDBFile

def fixed_pdb_file(input_pdb_file: Path or str, output_pdb_file= None)  -> Path:
    """
    Generate a fixed PDB file from the given PDB file.

    Parameters:
        pdb_file (Path or str): The path to the input PDB file.
        output_file (Path or str, optional): The path to the output PDB file. If not provided, the output file will have the same name as the input file.

    Returns:
        Path: The path to the output PDB file.

    """
    input_pdb_file = Path(input_pdb_file)
    fixer = PDBFixer(filename=str(input_pdb_file))

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()

    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)    
    
    if output_pdb_file is not None:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_file, 'w'))
    else:
        output_pdb_file=input_pdb_file
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_file, 'w'))
    return output_pdb_file

if __name__ == "__main__" :
    input_pdb_file = "/mnt/nas1/lanwei-125/MC5R/dock/hpep-complex/DEQPL/DEQPL.pdb"
    fixed_pdb_file(input_pdb_file)