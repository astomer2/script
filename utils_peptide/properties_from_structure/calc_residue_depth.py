from Bio.PDB import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth
from dataclasses import dataclass
import numpy as np
import csv
import os
from utils_comm.log_util import logger
from utils_peptide.properties_from_structure.identify_interface_residues import identify_interface_residues

MSMS_EXEC = "/mnt/nas/yuliu/softwares/msms/msms.x86_64Linux2.2.6.1"


@dataclass
class ResidueDepthInfo:
    """
    A class to represent the depth information of a residue.

    Attributes
    ----------
    pdb_file : str
        The name of the pdb file that the residue belongs to.
    chain_id : str
        The ID of the chain that the residue belongs to.
    residue_id : tuple
        The ID of the residue.
    residue_depth : float
        The depth of the residue.
    ca_depth : float
        The depth of the alpha carbon atom in the residue.
    """

    pdb_file: str
    chain_id: str
    residue_id: tuple
    residue_depth: float
    ca_depth: float

def interface_residues_depth(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Get the interacting residues for the current pdb_file
    interaction_result = interface_residues_results.get(pdb_path, None)

    if interaction_result is not None:
        interacting_residues = interaction_result.interacting_chainA_res
    else:
        interacting_residues = []

    for model in structure:
        # Calculate the residue depth
        rd = ResidueDepth(model, msms_exec=MSMS_EXEC)
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_id = residue.get_id()

                # Only calculate depth for interacting residues
                if np.isin(res_id, interacting_residues).any():
                    RD = rd[chain_id, res_id]

    return ResidueDepthInfo(pdb_path, chain_id, res_id, RD[0], RD[1])

def cal_residue_depth(workdir, interface_residues_results):
    """
    Calculate the depth of each residue in a protein structure.

    Parameters
    ----------
    workdir : str
        The directory that contains the pdb files of the protein structures to be analyzed.
    interface_residues_results : dict
        The output from the identify_interface_residues function. The dictionary's keys are PDB filenames, 
        and the values are the corresponding 'InteractionResult' objects, which contain information about 
        the interacting residues for that protein.

    Returns
    -------
    list of ResidueDepthInfo
        The depth information of each residue of the target protein in the interface.
    """
    total_pdb_files = [f for f in os.listdir(workdir) if f.endswith(".pdb")]
    total_pdb_count = len(total_pdb_files)
    logger.info(f"Total number of PDB files: {total_pdb_count}")

    # Initialize the PDB parser
    parser = PDBParser(QUIET=True)
    depths = []
    processed_pdb_count = 0
    # Loop over all pdb files in the work directory
    for pdb_file in os.listdir(workdir):
        # Skip files that are not pdb files
        if not pdb_file.endswith(".pdb"):
            continue
        processed_pdb_count += 1
        logger.info(
            f"Processing PDB file: {pdb_file} ({processed_pdb_count}/{total_pdb_count})"
        )
        pdb_path = os.path.join(workdir, pdb_file)
        # Parse the pdb file
        depths.append(interface_residues_depth(pdb_path))
    return depths




def write_results_to_csv(depths, csv_file):
    """
    Write the depth information of residues to a CSV file.

    Parameters
    ----------
    depths : list of ResidueDepthInfo
        The depth information of each residue.
    csv_file : str
        The path to the CSV file to be written.
    """
    with open(csv_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["pdb_file", "Chain ID", "Residue ID", "Residue Depth", "Ca Depth"]
        )
        for depth_info in depths:
            writer.writerow(
                [
                    depth_info.pdb_file,
                    depth_info.chain_id,
                    depth_info.residue_id,
                    depth_info.residue_depth,
                    depth_info.ca_depth,
                ]
            )


if __name__ == "__main__":
    workdir = (
        "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/cyc_pep_protein_complex"
    )
    logger.info("start the interface residues identification")
    interface_residues_results = identify_interface_residues(workdir)
    logger.info("finish identifying interface residues")

    logger.info("start the residue depth calculation")
    depths = cal_residue_depth(workdir, interface_residues_results)
    logger.info("finish calculating residue depths")

    csv_file = "/mnt/nas/yuliu/repos/peptide-deploy/test_outputs/residue_depths.csv"
    write_results_to_csv(depths, csv_file)
    logger.info("finish generating results in .csv")
