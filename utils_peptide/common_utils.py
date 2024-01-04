import os


def read_pdb_chain(workdir):
    """
    Read PDB files in the working directory and return the chains.

    Args:
    workdir (str): the working directory that contains the PDB files.

    Returns:
    dict: a dictionary where the keys are tuples of (pdb_name, chain_id) and the values are lists of lines in the PDB file.
    """
    chains = {}
    # for pdb_name in [f for f in os.listdir(workdir) if f.endswith(".pdb")]:
    for pdb_name in [
        f
        for f in os.listdir(workdir)
        if f.endswith(".pdb")
        and (not f.endswith("_A.pdb") and not f.endswith("_B.pdb"))
    ]:
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
