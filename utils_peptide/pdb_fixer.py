import os
from pdbfixer import PDBFixer
from pathlib import Path
from openmm.app import PDBFile

def fixed_pdb_file(pdb_file: Path or str, output_file= None)  -> Path:
    """
    Generate a fixed PDB file from the given PDB file.

    Parameters:
        pdb_file (Path or str): The path to the input PDB file.
        output_file (Path or str, optional): The path to the output PDB file. If not provided, the output file will have the same name as the input file.

    Returns:
        Path: The path to the output PDB file.

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
    
    if output_file is not None:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
    else:
        output_file=pdb_file
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
    return output_file

if __name__ == "__main__" :
    pass