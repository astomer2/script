import os
from pdbfixer import PDBFixer
from pathlib import Path
from openmm.app import PDBFile

def fixed_pdb_file(pdb_file) :
    """
    Generate a fixed PDB file.

    Args:
        pdb_file (str): The path to the PDB file.
        workdir (str): The working directory where the fixed PDB file will be saved.

    Returns:
        a fixed PDB file with the same name as the input PDB file.
    """
    pdb_file = Path(pdb_file)


    fixer = PDBFixer(filename= str(pdb_file))

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()

    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)    
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'))

    return pdb_file