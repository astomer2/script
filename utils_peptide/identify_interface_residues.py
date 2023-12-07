"""
This script identifies interface residues in protein-peptide complexes 
and writes the results to a CSV file when testing.
"""
import os
import csv

from prody.measure.contacts import findNeighbors
from prody import parsePDB
from utils_peptide.common_utils import read_pdb_chain
from utils_peptide.common_utils import write_pdb_chain
from dataclasses import dataclass


@dataclass
class InteractionResult:
    """
    A data class that stores the interaction result of a protein-peptide complex.

    Attributes:
        protein (str): The name of the protein file.
        peptide (str): The name of the peptide file.
        interacting_chainA_res (list): A list of interacting residues on chain A.
    """

    protein: str
    peptide: str
    interacting_chainA_res: list


def identify_interface_residues(workdir, needs_b=0):
    """
    This function identifies the interface residues in protein-peptide complexes.

    Args:
        workdir (str): The working directory that contains the PDB files.
        needs_b (int, optional): If set to 1, the function will also find the interacting residues on chain B. Defaults to 0.

    Returns:
        dict: A dictionary where the key is a tuple of (protein, peptide), and the value is an instance of InteractionResult.
    """

    result_dict = {}
    os.chdir(workdir)

    # Create a subdirectory named 'output'
    output_dir = os.path.join(workdir, "output")
    os.makedirs(output_dir, exist_ok=True)

    chains = read_pdb_chain(workdir)
    # Pass the output directory to write_pdb_chain
    pdb_filenames = write_pdb_chain(chains, output_dir)

    # Assuming the list has even number of elements
    for protein, peptide in zip(pdb_filenames[::2], pdb_filenames[1::2]):
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
        result_dict[(protein, peptide)] = InteractionResult(
            protein, peptide, interacting_chainA_res
        )

    return result_dict


def write_results_to_csv(results, filename):
    """
    This function writes the interaction results to a CSV file.

    Args:
        results (dict): The interaction results, which is a dictionary returned by identify_interface_residues().
        filename (str): The output CSV filename.
    """

    with open(filename, "w", newline="") as csvfile:
        fieldnames = ["protein", "peptide", "interacting_chainA_res"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results.values():
            writer.writerow(result.__dict__)


if __name__ == "__main__":
    workdir = (
        "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/cyc_pep_protein_complex"
    )
    results = identify_interface_residues(workdir)
    write_results_to_csv(
        results,
        "/mnt/nas/yuliu/repos/peptide-deploy/test_outputs/identify_interface_results.csv",
    )
