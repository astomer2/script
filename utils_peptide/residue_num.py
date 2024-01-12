def get_residue_num(pdb_path):
    """
    Counts the number of unique residues in a PDB file.

    Args:
        pdb_path (str): The path to the PDB file.

    Returns:
        int: The number of unique residues in the PDB file.
    """
    with open(pdb_path) as f:
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
    return residue_count

if __name__ == "__main__":
    import sys
    pdb_path = sys.argv[1]
    print(get_residue_num(pdb_path))