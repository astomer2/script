import os
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from dataclasses import dataclass
from typing import List
import csv


def read_pdb_chain(workdir):
    """
    Read PDB files in the working directory and return the chains.

    Args:
    workdir (str): the working directory that contains the PDB files.

    Returns:
    dict: a dictionary where the keys are tuples of (pdb_name, chain_id) and the values are lists of lines in the PDB file.
    """
    chains = {}
    for pdb_name in [f for f in os.listdir(workdir) if f.endswith(".pdb")]:
        # Remove the '.pdb' suffix
        pdb_base_name = pdb_name.split(".")[0]
        input_pdb = os.path.join(workdir, pdb_name)

        with open(input_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    chains.setdefault((pdb_base_name, chain_id), []).append(line)
    return chains


def write_pdb_chain(chains):
    """
    Write chains to separate PDB files.

    Args:
    chains (dict): a dictionary where the keys are tuples of (pdb_name, chain_id) and the values are lists of lines in the PDB file.
    
    Returns:
    list: a list of filenames of the written PDB files.
    """
    filenames = []
    for (pdb_name, chain_id), lines in chains.items():
        filename = f"{pdb_name}_{chain_id}.pdb"
        with open(filename, "w") as f:
            for line in lines:
                f.write(line)
            f.write("TER\n")
        filenames.append(filename)
    return filenames



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
    chains = read_pdb_chain(workdir)
    write_pdb_chain(chains)

    processed_pdb_names = set()
    results = []

    for (pdb_name, chain_id), _ in chains.items():
        # Skip if pdb_name has already been processed to avoid processing the same file for two different chains.
        if pdb_name in processed_pdb_names:
            continue

        protein_filename = f"{pdb_name}_A.pdb"
        peptide_filename = f"{pdb_name}_B.pdb"
        complex_filename = f"{pdb_name}.pdb"

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
    workdir = "/mnt/nas1/lanwei-125/FGF5/test"
    results = calc_relative_SASA(workdir)
    write_results_to_csv(results, "/mnt/nas1/lanwei-125/FGF5/test/SASA_results.csv")