import os
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO
from Bio.SeqUtils import seq1
from utils_comm.log_util import ic, logger


# https://www.wikidoc.org/index.php/Secondary_structure#The_DSSP_code
# G = 3-turn helix (310 helix). Min length 3 residues.
# H = 4-turn helix (alpha helix). Min length 4 residues.
# I = 5-turn helix (pi helix). Min length 5 residues.
other_helices_codes = ('G', 'I')


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


def write_pdb_chain(chains, sub_dir):
    """
    Write chains to separate PDB files in the given subdirectory.

    Args:
    chains (dict): a dictionary where the keys are tuples of (pdb_name, chain_id) and the values are lists of lines in the PDB file.
    sub_dir (str): the subdirectory where the PDB files will be written.

    Returns:
    list: a list of filenames of the written PDB files.
    """
    filenames = []
    for (pdb_name, chain_id), lines in chains.items():
        filename = f"{pdb_name}_{chain_id}.pdb"
        # Write the PDB file to the subdirectory
        output_path = os.path.join(sub_dir, filename)
        with open(output_path, "w") as f:
            for line in lines:
                f.write(line)
            f.write("TER\n")
        filenames.append(output_path)
    return filenames


def get_secondary_structure_from_pdb_file(pdb_file):
    """  """
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    # Calculate the secondary structure using DSSP
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='/mnt/nas1/bio-software/dssp-2.0.4-linux-amd64')
    return dssp


def get_helix_ratio_from_pdb_file(pdb_file):
    dssp = get_secondary_structure_from_pdb_file(pdb_file)
    alpha_helices_aa_count = 0
    other_helices_aa_count = 0
    all_aa_count = 0
    logger.info('%s', pdb_file)
    for residue in dssp:
        logger.info(residue[:3])
        if residue[2] == 'H':
            alpha_helices_aa_count += 1
        elif residue[2] in other_helices_codes:
            other_helices_aa_count += 1
        all_aa_count += 1
    assert all_aa_count > 0
    alpha_helices_aa_ratio = alpha_helices_aa_count / all_aa_count
    other_helices_aa_ratio = other_helices_aa_count / all_aa_count
    helices_aa_ratio = (alpha_helices_aa_count + other_helices_aa_count) / all_aa_count
    results = {
        'alpha_helices_aa_ratio': alpha_helices_aa_ratio,
        'other_helices_aa_ratio': other_helices_aa_ratio,
        'helices_aa_ratio': helices_aa_ratio,
        'all_aa_count': all_aa_count
    }
    return results


def check_alpha_helix_peptide_pdb(pdb_out_dir: Path, pep_chain="B"):
    """ Assert the peptide chain is "B". 
    """
    pep_dir = pdb_out_dir / "pep"
    pep_dir.mkdir(parents=True, exist_ok=True)
    pep_with_helix = set()
    for pdb_file in pdb_out_dir.glob("*.pdb"):
        pep_pdb_file = pep_dir / pdb_file.name
        if 'helix-ratio-' in pdb_file.stem:
            pep_with_helix.add(pdb_file.name)
            continue
        split_peptide_chain_from_pdb(pdb_file, pep_pdb_file, pep_chain)
        helix_ratio_result = get_helix_ratio_from_pdb_file(pep_pdb_file)
        helices_aa_ratio = helix_ratio_result["helices_aa_ratio"]
        alpha_helices_aa_ratio = helix_ratio_result["alpha_helices_aa_ratio"]
        if helices_aa_ratio > 0:
            if alpha_helices_aa_ratio > 0:
                helix_type = 'alpha_helix'
            else:
                helix_type = 'other_helix'
            new_pdb_file = pdb_out_dir / f"{pdb_file.stem}-{helix_type}-ratio-{helices_aa_ratio:.2f}.pdb"
            pep_with_helix.add(new_pdb_file.name)
            if not new_pdb_file.is_file():
                os.rename(pdb_file, new_pdb_file)
            # pdb_file and new_pdb_file is the same file, so delete the pdf_file.
            elif pdb_file.is_file():
                pdb_file.unlink()
    pep_with_helix = sorted(pep_with_helix)
    ic(pep_with_helix, len(pep_with_helix))


def split_peptide_chain_from_pdb(complex_pdb_file, sub_chain_pdb_file, saved_chain="B"):
    pdb_parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = pdb_parser.get_structure("reference", complex_pdb_file)
    model = structure[0]
    with open(sub_chain_pdb_file, "w", encoding='utf-8') as f:
        for chain in model:
            if chain.id == saved_chain:
                io.set_structure(chain)
                io.save(f)
                break


def calculate_protein_structural_similarity_by_rmsd(protein1 = 'protein1.pdb', protein2 = 'protein2.pdb'):
    """ tm-score is better than rmsd """
    # Load the two protein structures
    parser = PDBParser()
    structure1 = parser.get_structure('protein1', protein1)
    structure2 = parser.get_structure('protein2', protein2)

    # Align the two structures
    aligner = Superimposer()
    aligner.set_atoms(structure1.get_atoms(), structure2.get_atoms())
    aligner.apply(structure2.get_atoms())

    # Calculate the RMSD
    rmsd = aligner.rms
    logger.info(f"RMSD: {rmsd:.2f} Å")
    return rmsd


def create_pdb_from_partial_aas_from_pdb(residue_start_id: int, residue_end_id: int, pdb_file, output_pdb_file):
    """ Including residue_end_id """
    # Load the protein structure
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Select the amino acids with residue IDs between 10 and 20
    selected_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue_start_id <= residue.id[1] <= residue_end_id:
                    selected_residues.append(residue)

    # Create a new structure with only the selected amino acids
    selected_structure = structure.copy()
    for model in selected_structure:
        for chain in model:
            for residue in chain:
                if residue not in selected_residues:
                    chain.detach_child(residue.id)

    # Write the new structure to a PDB file
    io = PDBIO()
    io.set_structure(selected_structure)
    io.save(output_pdb_file)


def check_pdb_chain(pdb_file):
    """ check pdb by chain """
    for record in SeqIO.parse(pdb_file, "pdb-seqres"):
        pdb_id = record.id.strip().split(':')[0]
        chain  = record.annotations["chain"]
        ic(pdb_id, chain, record.seq[:10])
        if len(record.dbxrefs):
            _, UNP_id = record.dbxrefs[0].strip().split(':')
            ic(UNP_id)

def get_sequence_from_chain(PDB_file_path:str or Path,target_chain_id='A') -> str:
    '''get chain sequence with one letter code'''
    PDB_file_path = Path(PDB_file_path)
    pdbparser = PDBParser()
    structure = pdbparser.get_structure("X", PDB_file_path)
    chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
    query_chain_seq = chains[target_chain_id]
    return query_chain_seq

if __name__ == "__main__":
    file = ''
    get_secondary_structure_from_pdb_file(file)
