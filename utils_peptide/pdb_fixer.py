import os
from pdbfixer import PDBFixer
from pathlib import Path
from simtk.openmm.app import PDBFile

def fixed_pdb_file(pdb_file, workdir) :
    """
    Generate a fixed PDB file.

    Args:
        pdb_file (str): The path to the PDB file.
        workdir (str): The working directory where the fixed PDB file will be saved.

    Returns:
        a fixed PDB file with the same name as the input PDB file.
    """
    complex_fixed_pdb = Path(workdir) / Path(pdb_file).name


    fixer = PDBFixer(filename= str(complex_fixed_pdb))

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()

    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    PDBFile.writeFile(fixer.topology, fixer.positions, open(complex_fixed_pdb, 'w'))

    return complex_fixed_pdb