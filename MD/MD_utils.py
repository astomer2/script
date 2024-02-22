from itertools import cycle
import os
from pathlib import Path
from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO

def extract_three_letter_code(peptide_structure):
    """
    Extracts the three-letter amino acid sequences, sequence lengths, and file names from .pdb files in a given directory.

    :param work_path: Path to the directory containing the .pdb files.
    :return: A tuple of three lists: sequences, sequence_lengths, and file_names.
    """
    sequences = []

    with open(peptide_structure, "r", encoding="utf-8") as f:
        lines = f.readlines()
    sequence = []
    residue_ids = []
    for line in lines:
        if line.startswith("ATOM") :
            res_id = line[22:27].strip()
            if res_id not in residue_ids:
                residue_ids.append(res_id)
                res_name = line[17:20]
                sequence.append(res_name)
    sequence_str = " ".join(sequence)
    sequences.append(sequence_str)
    sequence_lengths = len(sequence)
    return sequences, sequence_lengths


def reference_id(protein_structure):
    """
    Generates a reference ID based on the protein structure.
    Reads the protein structure from the given file path and retrieves the lines.
    Iterates through the lines and extracts the protein ID from each line that starts
    with 'ATOM'. If the protein ID is different from the previous one, increments the
    residue count. Finally, constructs the reference ID string and returns the protein ID
    and reference ID.
    Parameters:
        protein_structure (str): The file path of the protein structure.
    Returns:
        tuple: A tuple containing the protein ID (str) and reference ID (str).
    """
    with open(protein_structure) as f:
        lines = f.readlines()
    protein_id = 0
    residue_count = 0
    prev_res_id = 0
    for line in lines:
        if line.startswith('ATOM'):
            protein_id = line[22:26].strip()
            if protein_id != prev_res_id:
                residue_count += 1
            prev_res_id = protein_id
    reference = f"1-{protein_id}"
    print(reference)
    print(protein_id)
    return protein_id, reference

def target_id(protein_id, peptide_structure):
    """
    Function to generate a target ID based on a protein ID and peptide structure.

    Parameters:
    protein_id (int): The ID of the protein.

    Returns:
    str: The generated target ID.

    This function reads the peptide structure from a file and calculates the target ID based on the protein ID and peptide ID. It iterates over the lines in the file, extracting the peptide ID from lines starting with 'ATOM'. If the peptide ID is different from the previous one, the residue count is incremented. Finally, the target ID is calculated as the sum of the protein ID and the peptide ID, plus 1.

    Example:
    >>> target_id(5)
    '6-11'
    """
    with open(peptide_structure) as f:
        lines = f.readlines()

    residue_count = 0
    prev_res_id = 0
    peptide_id = 0 
    for line in lines:
        if line.startswith('ATOM'):
            peptide_id = line[22:26].strip()
            if peptide_id != prev_res_id:
                residue_count += 1
            prev_res_id = peptide_id

    target = f"{int(protein_id)+1}-{int(peptide_id)+int(protein_id)}"
    print(peptide_id)
    print (target)
    return target

def identify_cys_id(peptide_structure):
    """
    This function takes in the path of a peptide structure file and identifies the residue IDs of all CYS or CYX residues in the file.

    Parameters:
    - peptide_structure (str): The path of the peptide structure file.

    Returns:
    - residue_ids (list): A list of integer residue IDs corresponding to the CYS or CYX residues found in the file.
    """
    cys_residue_ids = []
    prev_id = 0
    with open (peptide_structure) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('ATOM'):
            residue = line[17:20]
            if residue == 'CYS' or residue == 'CYX':
                residue_id = line[22:26].strip()
                if residue_id != prev_id:
                    prev_id = residue_id
                    cys_residue_ids.append(int(residue_id))
    return cys_residue_ids

def change_cys_to_cyx(peptide_structure):
    """
    Change all occurrences of 'CYS' to 'CYX' in the given peptide structure file.

    Parameters:
    - peptide_structure (str): The path to the peptide structure file.

    Returns:
    - None

    Raises:
    - FileNotFoundError: If the peptide structure file does not exist.
    """
    new_line = []
    with open (peptide_structure) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('ATOM'):
            residue = line[17:20]
            if residue == 'CYS' or residue == 'CYX':
                line = line[:17] + 'CYX' + line[20:]
                new_line.append(line)
            else :
                new_line.append(line)
    with open(peptide_structure, 'w') as f:
        f.writelines(new_line)
    return 


def check_last_residue_for_OXT(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]  # Assuming there's only one model in the structure
    last_residue = list(model.get_residues())[-1]
    has_OXT_atom = False
    for atom in last_residue:
        print(atom.get_name())
        if atom.get_name() == "OXT":
            has_OXT_atom = True
            break

    return has_OXT_atom
