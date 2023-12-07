from Bio.PDB import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth
from dataclasses import dataclass
import csv

MSMS_EXEC = "/mnt/nas/yuliu/softwares/msms/msms.x86_64Linux2.2.6.1"


@dataclass
class ResidueDepthInfo:
    """
    A class to represent the depth information of a residue.

    Attributes
    ----------
    chain_id : str
        The ID of the chain that the residue belongs to.
    residue_id : tuple
        The ID of the residue.
    residue_depth : float
        The depth of the residue.
    ca_depth : float
        The depth of the alpha carbon atom in the residue.
    """

    chain_id: str
    residue_id: tuple
    residue_depth: float
    ca_depth: float


def calculate_residue_depth(pdb_file):
    """
    Calculate the depth of each residue in a protein structure.

    Parameters
    ----------
    pdb_file : str
        The path to the PDB file of the protein structure.
    msms_exec : str
        The path to the executable of MSMS.

    Returns
    -------
    list of ResidueDepthInfo
        The depth information of each residue.
    """

    # Parse PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    depths = []
    for model in structure:
        rd = ResidueDepth(model, msms_exec=MSMS_EXEC)
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_id = residue.get_id()

                RD1 = rd[chain_id, res_id]
                # print("Chain ID : %s, Residue ID : %s, Residue depth : %s, Ca depth : %s" % (chain_id, res_id, RD1[0], RD1[1]))
                depths.append(ResidueDepthInfo(chain_id, res_id, RD1[0], RD1[1]))

    return depths


def write_depths_to_csv(depths, csv_file):
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
        writer.writerow(["Chain ID", "Residue ID", "Residue Depth", "Ca Depth"])
        for depth_info in depths:
            writer.writerow(
                [
                    depth_info.chain_id,
                    depth_info.residue_id,
                    depth_info.residue_depth,
                    depth_info.ca_depth,
                ]
            )


if __name__ == "__main__":
    pdb_file = "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/cyc_pep_protein_complex/test.pdb"
    depths = calculate_residue_depth(pdb_file)
    write_depths_to_csv(
        depths, "/mnt/nas/yuliu/repos/peptide-deploy/test_outputs/depths.csv"
    )
