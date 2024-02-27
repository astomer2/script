import logging
import time
from pathlib import Path
from Bio.PDB import PDBParser, Superimposer, PDBIO
import argparse

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def read_reference_pdb( reference_chain_id, reference_path):
    """
    Reads a reference PDB file and extracts the specified atoms to be aligned.

    Parameters:
    - reference_chain_id (str): The chain ID of the reference chain.
    - reference_path (str): The file path of the reference PDB file.
    - atoms_to_be_aligned (list): A list of atom IDs to be aligned.

    Returns:
    - ref_atoms (list): A list of reference atoms to be aligned.
    - ref_model (Model): The reference model containing the atoms.

    """
    reference_path = Path(reference_path)
    pdb_parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = pdb_parser.get_structure("reference", reference_path)
    model = structure[0]

    ref_atoms = []    
    for chain in model:
        if chain.get_id() ==  reference_chain_id:
            for res in chain:
                ref_atoms.append(res['CA'])
    ref_receptor = reference_path.parent/'reference.pdb'      
    with open(ref_receptor, 'w') as f:
        for ref_chains in model:
            if ref_chains.id == reference_chain_id:
                io.set_structure(ref_chains)
                io.save(f)
    return ref_atoms

def read_alignment_pdb( reference_chain_id, alignment_pdb):
    """
    Read an alignment PDB file and extract the atoms to be aligned.

    Parameters:
        reference_chain_id (str): The chain ID of the reference chain.
        alignment_pdb (str): The path to the alignment PDB file.
        atoms_to_be_aligned (list): A list of atom IDs to be aligned.

    Returns:
        tuple: A tuple containing the alignment atoms and the alignment model.
            - alignment_atoms (list): A list of alignment atoms.
            - alignment_model (Model): The alignment model.
    """

    pdb_parser = PDBParser(QUIET=True)
    alignment_structure = pdb_parser.get_structure("sample", alignment_pdb)
    alignment_model = alignment_structure[0]
    alignment_atoms = []    
    for chain in alignment_model:
        if chain.get_id() ==  reference_chain_id:
            for res in chain:
                alignment_atoms.append(res['CA'])
    return alignment_atoms, alignment_model

def superimpose_atoms(chain_id, ref_atoms, alignment_atoms, alignment_model, output_folder, alignment_pdb):
    """
    Superimposes atoms between two protein structures and saves the aligned structure as a PDB file.

    Parameters:
    - chain_id (str): The chain ID of the alignment chains.
    - ref_atoms (list): The reference atoms for superimposition.
    - alignment_atoms (list): The alignment atoms for superimposition.
    - alignment_model (object): The alignment model containing the atoms.
    - output_folder (str): The folder to save the aligned structure.
    - alignment_pdb (str): The name of the alignment PDB file.

    Returns:
    None
    """

    io = PDBIO()
    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, alignment_atoms)
    superimposer.apply(alignment_model.get_atoms())

    alignment_pdb_path = output_folder / f"{alignment_pdb.stem}{alignment_pdb.suffix}"

    with open(alignment_pdb_path, 'w') as fhandle:
        for alignment_chains in alignment_model:
            if alignment_chains.id == chain_id:
                io.set_structure(alignment_chains)
                io.save(fhandle)

    now_times = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    logger.info(f"{now_times} Superimpose alignment pdb file saved to {alignment_pdb_path}")

def alignment_pdb_from_reference( reference_chain_id, reference_path, pdb_path, aligment_chain_id, output_folder):
    ref_atoms = read_reference_pdb( reference_chain_id, reference_path)

    # Check if the output folder exists, create if not
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    for alignment_pdb in Path(pdb_path).glob("*.pdb"):
        if alignment_pdb.stem != 'reference':
            alignment_atoms, alignment_model = read_alignment_pdb( reference_chain_id, alignment_pdb)
            superimpose_atoms(aligment_chain_id, ref_atoms, alignment_atoms, alignment_model, output_folder, alignment_pdb)

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description='Multi-structure alignment,please ensure that the reference structure entered is the same as the receptor for which the structure needs to be aligned.')
    parse.add_argument('-ref','--reference_path', type=str,help='input the path of you selected reference structure')
    parse.add_argument('-i','--input_pdb_path', type=str, help='input the dir path of your alignment structures')
    parse.add_argument('-o','--output_folder', type=str,  help='the path of your output folder,the output files is peptide not complex')
    parse.add_argument('-ali_id','--alignment_chain_id', type=str,  help='the chain id of your alignment structure')
    parse.add_argument('-ref_id','--reference_chain_id', type=str,  help='the chain id of your reference structure')
    
    args = parse.parse_args()
    reference_path = args.reference_path
    pdb_path = args.input_pdb_path
    output_folder = args.output_folder
    aligment_chain_id = args.alignment_chain_id
    reference_chain_id = args.reference_chain_id

    alignment_pdb_from_reference(reference_chain_id, reference_path, pdb_path, aligment_chain_id, output_folder)