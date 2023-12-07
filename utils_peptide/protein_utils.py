from Bio.PDB import PDBParser, DSSP, Superimposer, PDBIO
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=str)
ic.lineWrapWidth = 120


def get_secondary_structure_from_pdb_file(pdb_file):
    """  """
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Calculate the secondary structure using DSSP
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='/mnt/nas1/bio-software/dssp-2.0.4-linux-amd64')

    # Print the secondary structure assignments for each residue
    for residue in dssp:
        print(residue)


def calculate_protein_structural_similarity():
    """ fake code, tm-score is better than rmsd """
    # Load the two protein structures
    parser = PDBParser()
    structure1 = parser.get_structure('protein1', 'protein1.pdb')
    structure2 = parser.get_structure('protein2', 'protein2.pdb')

    # Align the two structures
    aligner = Superimposer()
    aligner.set_atoms(structure1.get_atoms(), structure2.get_atoms())
    aligner.apply(structure2.get_atoms())

    # Calculate the RMSD
    rmsd = aligner.rms
    print(f"RMSD: {rmsd:.2f} Å")
    return rmsd


def create_pdb_from_partial_aas_from_pdb():
    """ fake code """
    # Load the protein structure
    parser = PDBParser()
    structure = parser.get_structure('protein', 'protein.pdb')

    # Select the amino acids with residue IDs between 10 and 20
    selected_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 10 <= residue.id[1] <= 20:
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
    io.save('selected_amino_acids.pdb')


if __name__ == "__main__":
    file = '/mnt/sdd/alphafold/af_out/peptides/DTFGRCRRWWAALGACRR/ranked_1.pdb'
    get_secondary_structure_from_pdb_file(file)