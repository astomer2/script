import logging
import time
from pathlib import Path
from Bio.PDB import PDBParser, Superimposer, PDBIO

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def aligned_atoms(start_id, end_id):
    atoms_to_be_aligned = range(start_id, end_id + 1)
    return atoms_to_be_aligned

def read_reference_pdb( reference_chain_id, reference_path, atoms_to_be_aligned):
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

    pdb_parser = PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reference_path)
    ref_model = ref_structure[0]
    ref_atoms = []

    for ref_chain in ref_model:
        if ref_chain.id ==  reference_chain_id:
            for ref_res in ref_chain:
                if ref_res.get_id()[1] in atoms_to_be_aligned:
                    ref_atoms.append(ref_res['CA'])
    return ref_atoms, ref_model

def read_alignment_pdb( reference_chain_id, alignment_pdb, atoms_to_be_aligned):
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

    for sample_chain in alignment_model:
        if sample_chain.id ==  reference_chain_id:
            for sample_res in sample_chain:
                if sample_res.get_id()[1] in atoms_to_be_aligned:
                    alignment_atoms.append(sample_res['CA'])
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

def alignment_pdb_from_reference(start_id, end_id, reference_chain_id, reference_path, pdb_path, chain_id, output_folder):
    atoms_to_be_aligned = aligned_atoms(start_id, end_id)
    ref_atoms, _ = read_reference_pdb( reference_chain_id, reference_path, atoms_to_be_aligned)

    # Check if the output folder exists, create if not
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    for alignment_pdb in Path(pdb_path).glob("*.pdb"):
        if alignment_pdb != reference_path:
            alignment_atoms, alignment_model = read_alignment_pdb( reference_chain_id, alignment_pdb, atoms_to_be_aligned)
            superimpose_atoms(chain_id, ref_atoms, alignment_atoms, alignment_model, output_folder, alignment_pdb)

if __name__ == "__main__":
    reference_path ='/mnt/nas1/lanwei-125/FGF5/dock_prepare/FGF5.pdb'
    pdb_path = '/mnt/nas1/lanwei-125/FGF5/FGF5_all_cyc/'

    start_id = 1
    end_id = 131
    reference_chain_id = 'A'

    aligment_chain_id = 'B'
    output_folder = '/mnt/nas1/lanwei-125/FGF5/disulfide_peptide'
    alignment_pdb_from_reference(start_id, end_id, reference_chain_id, reference_path, pdb_path, aligment_chain_id, output_folder)