"""
This script identifies interface residues in protein-peptide complexes 
and writes the results to a CSV file when testing.
"""
import os
import csv

from prody.measure.contacts import findNeighbors
from prody import parsePDB
from utils_comm.log_util import logger
from utils_peptide.pdb_utils import read_pdb_chain
from utils_peptide.pdb_utils import write_pdb_chain
from dataclasses import dataclass


@dataclass
class InteractionResult:
    """
    A data class that stores the interaction result of a protein-peptide complex.

    Attributes:
        pdb_file (str): The base name of the PDB file.
        interacting_chainA_res (list): A list of interacting residues on chain A.
    """

    pdb_file: str
    interacting_chainA_res: list


def identify_interface_residues(workdir, needs_b=0):
    """
    This function identifies the interface residues in protein-peptide complexes.

    Args:
        workdir (str): The working directory that contains the PDB files.
        needs_b (int, optional): If set to 1, the function will also find the interacting residues on chain B. Defaults to 0.

    Returns:
        dict: A dictionary where the key is the name of the PDB file, and the value is an instance of InteractionResult.
    """

    result_dict = {}
    os.chdir(workdir)

    # Create a subdirectory named 'complex_split_out_pdb'
    output_dir = os.path.join(workdir, "complex_split_out_pdb")
    os.makedirs(output_dir, exist_ok=True)
    # Count the total number of all pdb files under the directory folder
    total_pdb_files = [f for f in os.listdir(workdir) if f.endswith(".pdb")]
    total_pdb_count = len(total_pdb_files)
    logger.info(f"Total number of PDB files: {total_pdb_count}")

    chains = read_pdb_chain(workdir)
    # Pass the output directory to write_pdb_chain
    pdb_filenames = write_pdb_chain(chains, output_dir)
    processed_pdb_count = 0
    # Assuming the list has even number of elements
    for protein, peptide in zip(pdb_filenames[::2], pdb_filenames[1::2]):
        processed_pdb_count += 1
        pdb_name = os.path.basename(protein).split("_")[
            0
        ]  # Extract the pdb_name from the protein file name
        pdb_file = f"{pdb_name}.pdb"
        logger.info(
            f"Processing PDB file: {pdb_file} ({processed_pdb_count}/{total_pdb_count})"
        )
        rec = parsePDB(protein)
        pep = parsePDB(peptide)

        near_n = findNeighbors(rec, 5, pep)

        interacting_chainA_res = []
        if needs_b == 1:
            interacting_chainB_res = []
        for a1, a2, d in near_n:
            interacting_chainA_res.append(a1.getResindex())
            if needs_b == 1:
                interacting_chainB_res.append(a2.getResindex())

        interacting_chainA_res = list(set(interacting_chainA_res))

        interacting_chainA_res.sort()

        # Store the result in the dictionary
        result_dict[(pdb_file)] = InteractionResult(pdb_file, interacting_chainA_res)

    return result_dict


def write_results_to_csv(results, filename):
    """
    This function writes the interaction results to a CSV file.

    Args:
        results (dict): The interaction results, which is a dictionary returned by identify_interface_residues().
        filename (str): The output CSV filename.
    """

    with open(filename, "w", newline="") as csvfile:
        fieldnames = ["pdb_file", "interacting_chainA_res"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results.values():
            writer.writerow(result.__dict__)


if __name__ == "__main__":
    logger.info("start identifying the interface residues")
    workdir = (
        "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/cyc_pep_protein_complex"
    )
    results = identify_interface_residues(workdir)
    write_results_to_csv(
        results,
        "/mnt/nas/yuliu/repos/peptide-deploy/test_outputs/identify_interface_results.csv",
    )
    logger.info("end")
