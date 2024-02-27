import os
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from dataclasses import dataclass
from typing import List
from utils_comm.log_util import logger
from utils_peptide.pdb_utils import read_pdb_chain
from utils_peptide.pdb_utils import write_pdb_chain
import csv

@dataclass
class SASAResult:
    """
    A data class that represents the SASA result.

    Attributes:
    pdb_name (str): the name of the PDB file.
    protein_sasa (float): the SASA of the protein.
    peptide_sasa (float): the SASA of the peptide.
    complex_sasa (float): the SASA of the complex.
    relative_sasa (float): the relative SASA of the complex.
    """

    pdb_name: str
    protein_sasa: float
    peptide_sasa: float
    complex_sasa: float
    relative_sasa: float


def calc_SASA(PDB_file):
    """
    Calculate the Solvent Accessible Surface Area (SASA) for the given PDB file.

    Args:
    PDB_file (str): the path to the PDB file.

    Returns:
    float: the calculated SASA.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", PDB_file)
    sr = ShrakeRupley()
    sr.compute(structure, level="S")
    SASA = 0
    for residue in structure.get_residues():
        for atom in residue:
            SASA += atom.sasa
    return SASA


def write_results_to_csv(results: List[SASAResult], filename: str):
    """
    Write the given list of SASAResult objects to a CSV file.

    Args:
    results (List[SASAResult]): a list of SASAResult objects.
    filename (str): the output CSV file name.
    """
    with open(filename, "w", newline="") as csvfile:
        fieldnames = [
            "pdb_name",
            "protein_sasa",
            "peptide_sasa",
            "complex_sasa",
            "relative_sasa",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def calc_relative_SASA(workdir) -> List[SASAResult]:
    """
    Calculate the relative SASA for each PDB file in the given directory and write the results to a CSV file.

    Args:
    workdir (str): the working directory that contains the PDB files.

    Returns:
    List[SASAResult]: a list of SASAResult objects, each representing the SASA calculation result for one PDB file.

    Note:
    This function will change the current working directory to `workdir`.
    """
    os.chdir(workdir)

    # Create a subdirectory named 'complex_split_out_pdb'
    output_dir = os.path.join(workdir, "complex_split_out_pdb")
    os.makedirs(output_dir, exist_ok=True)
    total_pdb_files = [f for f in os.listdir(workdir) if f.endswith(".pdb")]
    total_pdb_count = len(total_pdb_files)
    logger.info(f"Total number of PDB files: {total_pdb_count}")
    chains = read_pdb_chain(workdir)

    # Pass the output directory to write_pdb_chain
    write_pdb_chain(chains, output_dir)
    processed_pdb_names = set()
    results = []
    processed_pdb_count = 0
    for (pdb_name, chain_id), _ in chains.items():
        # Skip if pdb_name has already been processed to avoid processing the same file for two different chains.
        if pdb_name in processed_pdb_names:
            continue
        processed_pdb_count += 1
        logger.info(f"Processing PDB file: {pdb_name} ({processed_pdb_count}/{total_pdb_count})")
        # Create the full file path for the protein PDB file in the output directory
        protein_filename = os.path.join(output_dir, f"{pdb_name}_A.pdb")

        # Create the full file path for the peptide PDB file in the output directory
        peptide_filename = os.path.join(output_dir, f"{pdb_name}_B.pdb")

        # Create the full file path for the complex PDB file in the working directory
        complex_filename = os.path.join(workdir, f"{pdb_name}.pdb")

        protein_sasa = calc_SASA(protein_filename)
        peptide_sasa = calc_SASA(peptide_filename)
        complex_sasa = calc_SASA(complex_filename)
        relative_sasa = complex_sasa - (protein_sasa + peptide_sasa)

        # Create a new SASAResult instance to store the results for this pdb.
        result = SASAResult(
            pdb_name, protein_sasa, peptide_sasa, complex_sasa, relative_sasa
        )
        results.append(result)

        processed_pdb_names.add(pdb_name)
    return results


if __name__ == "__main__":
    logger.info("start calculating the relative SASA")
    workdir = (
        "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/cyc_pep_protein_complex"
    )
    results = calc_relative_SASA(workdir)
    write_results_to_csv(
        results, "/mnt/nas/yuliu/repos/peptide-deploy/test_outputs/SASA_results.csv"
    )
    logger.info("end")